import bisect
import numpy as np
from openpyxl import load_workbook
import pandas as pd
from scipy import interpolate
import warnings


def check_repeated_headings(input_file, sheet):
    """[summary]

    Args:
        input_file ([type]): [description]
        sheet ([type]): [description]

    Raises:
        ValueError: [description]
    """

    # Get the columns names that user can define (except the first four ones that are fixed). The variable columns is a tuple.
    columns = list(
        sheet.iter_rows(
            min_row=3,
            max_row=3,
            min_col=sheet.min_column,
            max_col=sheet.max_column,
            values_only=True,
        )
    )[0][4:]
    # Buil dictionay exploiting dict comprehension: each key as the numer of repetitions of the column as the corresponding value.
    dict_colum = {column: columns.count(column) for column in columns}
    # Raise error message
    if max(list(dict_colum.values())) > 1:
        raise ValueError(
            f"ERROR! Different objects of the same kind ({sheet['A1'].value}) can not have the same ID.\nUser defines the following:\n{dict_colum.items()}.\nPlease check the headers in sheet {sheet.title} of file {input_file}"
        )


# End fuction check_repeated_headings.


def check_headers(cond, path_input, path_operation, sheet_input, sheet_operation):
    """[summary]

    Args:
        cond ([type]): [description]
        path_input ([type]): [description]
        path_operation ([type]): [description]
        sheet_input ([type]): [description]
        sheet_operation ([type]): [description]

    Raises:
        SyntaxError: [description]
    """
    header_input = list(
        pd.read_excel(
            path_input, sheet_name=sheet_input.title, skiprows=2, header=0, index_col=0
        ).columns
    )
    header_operation = list(
        pd.read_excel(
            path_operation,
            sheet_name=sheet_operation.title,
            skiprows=2,
            header=0,
            index_col=0,
        ).columns
    )

    for ii in range(len(header_input)):
        if header_input[ii] != header_operation[ii]:
            raise SyntaxError(
                f"ERROR in method {cond.__init__.__name__} of class {cond.__class__.__name__}: object headers in {sheet_input} of file {cond.file_input['STRUCTURE_ELEMENTS']} and in {sheet_operation} of file {cond.file_input['OPERATION']} should be the same!\nCompare the third row, column {ii+1} of those sheets.\n"
            )


# End function check_headers


def check_object_number(object, path_1, path_2, sheet_1, sheet_2):
    """[summary]

    Raises:
        ValueError: [description]
    """
    dict_method = dict(
        Simulation="conductor_instance", Conductor="conductor_components_instance"
    )
    if int(sheet_1.cell(row=1, column=2).value) != int(
        sheet_2.cell(row=1, column=2).value
    ):
        raise ValueError(
            f"ERROR in class {object.__class__.__name__} method {dict_method[object.__class__.__name__]}: number of objects defined in file {path_1} sheet {sheet_1.title} and in file {path_2} sheet {sheet_2.title} must be the same.\nPlease compare cell B1 of those sheets.\n"
        )


# End function check_object_number.


def set_diagnostic(vv, **kwargs):
    """[summary]

    Args:
        vv ([type]): [description]

    Returns:
        [type]: [description]
    """

    # Keep only the positive values, negative values means that user do not want to save spatial distribution (or time evolutions) of variables at that times (or spatial coordinate). Positive values are sorted in ascendent order.
    vv = np.sort(vv[vv >= 0])
    # add to array vv the lower boundary (lb) and the upper boundary (ub) to save data by default also at these times (or spatial coordinates)
    vv = np.concatenate((vv, np.array([kwargs["lb"], kwargs["ub"]], dtype=float)))
    return np.unique(vv)


# End function set_diagnostic.


def load_auxiliary_files(file_path, sheetname):
    """Function that load the auxiliary input file as a data frame

    Args:
        file_path (_type_): _description_
        sheetname (_type_): _description_

    Returns:
        _type_: _description_
    """
    wb = load_workbook(file_path, data_only=True)
    sheet = wb[sheetname]
    return (
        pd.read_excel(file_path, sheet_name=sheetname, header=0, index_col=0),
        sheet.cell(1, 1).value,
    )


# End function load_auxiliary_files


def build_interpolator(df, interpolation_kind="linear"):
    """Function that builds the interpolator for the data loaded from auxiliary input files. If data are costant in time or in space 1D interpolate.interp1d is used, if values describe a dependence in both time and space interpolate.interp1d is used.

    Args:
        df (_type_): _description_
        interpolation_kind (str, optional): _description_. Defaults to "linear".

    Returns:
        _type_: _description_
    """
    # Convert the index of the dataframe to space points used in the interpolation function.
    strand_space_points = df.index.to_numpy()

    # Convert the columns of the dataframe to time points used in the interpolation function.
    strand_time_points = df.columns.to_numpy()

    known_val = df.to_numpy()

    if strand_time_points.size == 1 and strand_space_points.size > 1:
        # Values are constant in time but not in space.
        known_val = known_val.reshape(known_val.shape[0])
        return (
            interpolate.interp1d(
                strand_space_points,
                known_val,
                bounds_error=False,
                fill_value=known_val[-1],
                kind=interpolation_kind,
            ),
            "space_only",
        )

    elif strand_time_points.size > 1 and strand_space_points.size == 1:
        # Values are constant in space but not in time.
        known_val = known_val.reshape(known_val.shape[1])
        return (
            interpolate.interp1d(
                strand_time_points,
                known_val,
                bounds_error=False,
                fill_value=known_val[-1],
                kind=interpolation_kind,
            ),
            "time_only",
        )

    elif strand_time_points.size > 1 and strand_space_points.size > 1:
        return (
            interpolate.interp2d(
                strand_time_points,
                strand_space_points,
                known_val,
                kind=interpolation_kind,
            ),
            "space_and_time",
        )


# End function build_interpolator


def do_interpolation(interpolator, zcoord, time_step, kind):
    """Functin that makes the interpolation of the data loaded from auxiliary input files accorting to the kind flag.

    Args:
        interpolator (_type_): _description_
        zcoord (_type_): _description_
        time_step (_type_): _description_
        kind (_type_): _description_

    Returns:
        _type_: _description_
    """

    if kind == "space_only":
        # Values are constant in time but not in space.
        return interpolator(zcoord)

    elif kind == "time_only":
        # Values are constant in space but not in time.
        return interpolator(time_step)

    elif kind == "space_and_time":
        return interpolator(time_step, zcoord).reshape(zcoord.shape)


# End function do_interpolation


def with_read_csv(fname, col_name, delimiter=";"):
    """[summary]

    Args:
        fname ([type]): [description]
        col_name ([type]): [description]
        delimiter (str, optional): [description]. Defaults to ";".

    Returns:
        [type]: [description]
    """
    # The following works both for a single file with the spatial discretization for all the conductors and for a file with the spatial discretization for each conductor (i.e. one file with a single column for each conductor). The header of the column must be f"prefix_{conductor.name} [SI units]".
    # Load the column of the file (exension dat, csv, or tsv) with the spatial discretization of the conductor named conductor.name as a series; then convert to numpy dropping NaN values.
    vector = (
        pd.read_csv(fname, sep=delimiter, usecols=[col_name], squeeze=True)
        .dropna()
        .to_numpy()
    )
    return vector, vector.shape[0]

    # End function with_read_csv.


def with_read_excel(fname, col_name, sheet=0):
    """[summary]

    Args:
        fname ([type]): [description]
        col_name ([type]): [description]
        sheet (int, optional): [description]. Defaults to 0.

    Returns:
        [type]: [description]
    """
    # The following works both for a single file with the spatial discretization for all the conductors and for a file with the spatial discretization for each conductor (i.e. one file with a single column for each conductor). The header of the column must be f"x_{conductor.name} [m]".
    # N.B. sheet non specified (to be added as features!!)
    # Load the column of the file .xlsv with the spatial discretization of the conductor named conductor.name as a series; then convert to numpy dropping NaN values.
    vector = (
        pd.read_excel(fname, sheet_name=sheet, usecols=[col_name], squeeze=True)
        .dropna()
        .to_numpy()
    )
    return vector, vector.shape[0]

    # End function with_read_csv.


def get_from_xlsx(conductor, f_path, comp, flag_name, *INTIAL):

    """
    ##############################################################################
    #              Get_from_xlsx(file_path, comp)
    ##############################################################################
    #
    # Function that perform interpolation on data read from .xlsx files to
    # initialize SolidComponent python objects parameters like magnetic field,
    # magnetic field gradient, strain and external heat.
    #
    ##############################################################################
    # VARIABLE    I/O    TYPE              DESCRIPTION                      UNIT
    # --------------------------------------------------------------------------
    # f_path      I      string            path of the file to be read      -
    # comp        I      object            python object of
    #                                      class SolidComponent            -
    # value       O      np array float    vector of interpolation
    #                                      results                          many
    # flag        O      scalar integer    flag that specifies data units   -
    ##############################################################################
    # Invoked functions/methods: Read_interp_file, interpolation
    #
    ##############################################################################
    # N.B. no self parameters in input.
    # N.B. In case of bfield.xlsx and alphaB.xlsx if flag == 1 value unit are T
    # and T/m respectively, if flag == 2 value is per current value (i.e T/A
    # and T/m/A).
    # N.B. value has the same shape of xcoord and is an np array
    ##############################################################################
    #
    # Author D. Placido Polito 06/2020
    #
    ##############################################################################
    """

    # Function Read_interp_file starts here. (cdp)
    # Read data from .xlsx file
    # INTIAL[0]
    if len(INTIAL) == 0:  # interpolation in time and space
        [flag, tt, matrix, xx, sheet] = read_interp_file(
            f_path, comp, action="time_&_space"
        )
        # Perform loop on mesh nodes and perform interpolation in both time and \
        # space.
        value = interpolation(
            conductor,
            comp,
            matrix,
            tt,
            f_path,
            sheet,
            xx,
            action="time_&_space",
            Flag_name=flag_name,
        )
    else:  # interpolation in time
        [flag, tt, matrix, sheet] = read_interp_file(
            f_path, comp, INTIAL[0], action="time"
        )
        # Perform loop on mesh nodes and perform interpolation in time.
        value = interpolation(conductor, comp, matrix, tt, f_path, sheet, action="time")
    return [value, flag]


def interpolation(conductor, comp, MM, tvec, f_path, sheet, *xvec, **options):

    """
    ##############################################################################
    #              Interpolation(conductor, MM, tvec, xvec, xcoord)
    ##############################################################################
    #
    # Function that interpolate on MM @ t = time and x = xcoord (cdp)
    #
    ##############################################################################
    # VARIABLE    I/O    TYPE                  DESCRIPTION                  UNIT
    # -------------------------------------------------------------------------
    # MM          I      np 2 dimensional      data to be interpolated
    #                    array float           read from file .xlsx         T\W
    # tvec        I      np array float        time vector @ which data
    #                                          are available                s
    # xvec        I      np array float        space vector @ which data
    #                                          are available                m
    # time        I      scalar float          time @ which interpolation
    #                                          is carried out               s
    # xcoord*     I      np array float        spatial coordinate @ which
    #                                          interpolation is carried out m
    # yy          O      np array float        interpolation result
    #                                          vector                       many
    ##############################################################################
    #
    # Invoked functions/methods: none
    #
    ##############################################################################
    # N.B yy has the same shape of xcoord and is an np array
    # * xcoord is given by conductor.xcoord
    ##############################################################################
    # Author D. Placido Polito 06/2020
    #
    ##############################################################################
    """

    # Function interpolation starts here. (cdp)

    if options.get("action") == "time_&_space":
        # time-space interpolation (cdp, 07/2020)
        # Check on xvec and tvec values: first value must be >= 0 and they should \
        # be sorted from min to max value (cdp, 06/2020)
        xvec = np.array(xvec[0])
        if (xvec[0] >= 0 and all(np.diff(xvec) >= 0)) and (
            tvec[0] >= 0 and all(np.diff(tvec) >= 0)
        ):
            # Check that time @ which I perform interpolation is within tvec values, \
            # if not yy is given by spatial interpolation @ time = tvec[0] or \
            # time = tvec[-1] (cdp, 07/2020)
            yy = np.zeros(conductor.gird_features["N_nod"], dtype=float)
            lower_bound = np.min(
                np.nonzero(conductor.gird_features["xcoord"] >= xvec[0])
            )
            upper_bound = np.max(
                np.nonzero(conductor.gird_features["xcoord"] <= xvec[-1])
            )
            if options["Flag_name"] == "IQFUN":
                # Interpolation to get external flux (cdp, 11/2020)
                if comp.operations["IQFUN"] == -1:
                    # square wave in time and space (cdp, 11/2020)
                    if (
                        conductor.cond_time[-1] >= tvec[0]
                        and conductor.cond_time[-1] < tvec[-1]
                        and conductor.cond_time[-1] >= 0
                    ):
                        # search in time (cdp, 11/2020)
                        # array smart (cdp, 11/2020)
                        ii = np.max(np.nonzero(tvec <= conductor.cond_time[-1]))
                        # Square wave in time and space: at any times xcoord < lower_bound \
                        # and xcoord > upper_bound are equal to 0.0 (cdp, 10/2020)
                        if conductor.cond_time[-1] <= tvec[-2]:
                            yy[lower_bound : upper_bound + 1] = MM[0, ii]
                        else:
                            yy[lower_bound : upper_bound + 1] = MM[0, ii + 1]
                        # end if (cdp, 10/2020)
                    # end if conductor.cond_time[-1] (cdp, 11/2020)
                # end if comp.operations["IQFUN"]
            else:
                if (
                    conductor.cond_time[-1] >= tvec[0]
                    and conductor.cond_time[-1] < tvec[-1]
                    and conductor.cond_time[-1] >= 0
                ):
                    yy[0:lower_bound] = MM[0, ii] + (
                        conductor.cond_time[-1] - tvec[ii]
                    ) / (tvec[ii + 1] - tvec[ii]) * (MM[0, ii + 1] - MM[0, ii])
                    for kk in range(len(xvec) - 1):
                        lb = np.min(
                            np.nonzero(
                                conductor.gird_features["xcoord"] >= xvec[kk]
                            )
                        )
                        ub = np.max(
                            np.nonzero(
                                conductor.gird_features["xcoord"] <= xvec[kk + 1]
                            )
                        )
                        if (ii < tvec.shape[0] - 1) and (ub <= upper_bound):
                            fx = (
                                conductor.gird_features["xcoord"][lb : ub + 1]
                                - xvec[kk]
                            ) / (xvec[kk + 1] - xvec[kk])
                            # interpolation in time (cdp)
                            yy_t_j = MM[kk, ii] + (
                                conductor.cond_time[-1] - tvec[ii]
                            ) / (tvec[ii + 1] - tvec[ii]) * (
                                MM[kk, ii + 1] - MM[kk, ii]
                            )
                            yy_t_jplus1 = MM[kk + 1, ii] + (
                                conductor.cond_time[-1] - tvec[ii]
                            ) / (tvec[ii + 1] - tvec[ii]) * (
                                MM[kk + 1, ii + 1] - MM[kk + 1, ii]
                            )
                            # interpolation in space (cdp)
                            yy[lb : ub + 1] = yy_t_j + (yy_t_jplus1 - yy_t_j) * fx
                            if ub == upper_bound:
                                yy[ub] = MM[kk + 1, ii] + (
                                    conductor.cond_time[-1] - tvec[ii]
                                ) / (tvec[ii + 1] - tvec[ii]) * (
                                    MM[kk + 1, ii + 1] - MM[kk + 1, ii]
                                )
                        elif (ii == tvec.shape[0] - 1) and (ub <= upper_bound):
                            fx = (
                                conductor.gird_features["xcoord"][lb : ub + 1]
                                - xvec[kk]
                            ) / (xvec[kk + 1] - xvec[kk])
                            yy[lb : ub + 1] = (
                                MM[kk, ii] + (MM[kk + 1, ii] - MM[kk, ii]) * fx
                            )
                            if ub == upper_bound:
                                yy[ub] = MM[kk + 1, ii]
                    yy[upper_bound + 1 : len(yy)] = MM[kk + 1, ii] + (
                        conductor.cond_time[-1] - tvec[ii]
                    ) / (tvec[ii + 1] - tvec[ii]) * (
                        MM[kk + 1, ii + 1] - MM[kk + 1, ii]
                    )
                elif conductor.cond_time[-1] < tvec[0] and conductor.cond_time[-1] >= 0:
                    # time > 0 always (cdp, 07/2020)
                    yy[0:lower_bound] = MM[0, 0]
                    for kk in range(len(xvec) - 1):
                        lb = np.min(
                            np.nonzero(
                                conductor.gird_features["xcoord"] >= xvec[kk]
                            )
                        )
                        ub = np.max(
                            np.nonzero(
                                conductor.gird_features["xcoord"] <= xvec[kk + 1]
                            )
                        )
                        fx = (
                            conductor.gird_features["xcoord"][lb : ub + 1]
                            - xvec[kk]
                        ) / (xvec[kk + 1] - xvec[kk])
                        # spatial interpolation @ time = tvec[0] (cdp, 07/2020)
                        yy[lb : ub + 1] = MM[kk, 0] + (MM[kk + 1, 0] - MM[kk, 0]) * fx
                    yy[upper_bound + 1 : len(yy)] = MM[kk + 1, 0]
                elif conductor.cond_time[-1] >= tvec[-1]:
                    # time > 0 always (cdp, 07/2020)
                    yy[0:lower_bound] = MM[0, -1]
                    for kk in range(len(xvec) - 1):
                        lb = np.min(
                            np.nonzero(
                                conductor.gird_features["xcoord"] >= xvec[kk]
                            )
                        )
                        ub = np.max(
                            np.nonzero(
                                conductor.gird_features["xcoord"] <= xvec[kk + 1]
                            )
                        )
                        fx = (
                            conductor.gird_features["xcoord"][lb : ub + 1]
                            - xvec[kk]
                        ) / (xvec[kk + 1] - xvec[kk])
                        # spatial interpolation @ time = tvec[-1] (cdp, 07/2020)
                        yy[lb : ub + 1] = (
                            MM[kk, -1] + (MM[kk + 1, -1] - MM[kk, -1]) * fx
                        )
                    yy[upper_bound + 1 : len(yy)] = MM[kk + 1, -1]
                else:
                    raise ValueError(
                        f"ERROR in {interpolation.__name__}: time must \
                        be positive scalar!\ntime = {conductor.cond_time[-1]}"
                    )
        else:
            raise ValueError(
                f"ERROR in {interpolation.__name__}: xvec and/or \
                         tvec are not sorted or have values < 0.0!\nCheck \
                         sheet {sheet} in {f_path}\n"
            )
    elif options.get("action") == "time":
        # time-only interpolation (cdp, 07/2020)
        if tvec[0] >= 0 and all(np.diff(tvec) >= 0):
            yy = np.zeros(len(MM[:, 0]))  # one value for each row on MM
            if (
                conductor.cond_time[-1] >= tvec[0]
                and conductor.cond_time[-1] < tvec[-1]
                and conductor.cond_time[-1] >= 0
            ):
                # search in time (cdp)
                # array smart (cdp)
                ii = np.max(np.nonzero(tvec <= conductor.cond_time[-1]))
                for kk in range(len(yy)):
                    yy[kk] = MM[kk, ii] + (conductor.cond_time[-1] - tvec[ii]) / (
                        tvec[ii + 1] - tvec[ii]
                    ) * (MM[kk, ii + 1] - MM[kk, ii])
            elif (
                conductor.cond_time[-1] < tvec[0] and conductor.cond_time[-1] >= 0
            ):  # time > 0 always (cdp, 07/2020)
                for kk in range(len(yy)):
                    yy[kk] = MM[kk, 0]
            elif conductor.cond_time[-1] >= tvec[-1]:
                # time > 0 always (cdp, 07/2020)
                for kk in range(len(yy)):
                    yy[kk] = MM[kk, -1]
            else:
                raise ValueError(
                    f"ERROR in {interpolation.__name__}: time must \
                        be positive scalar!\ntime = {conductor.cond_time[-1]}"
                )
        else:
            raise ValueError(
                f"ERROR in {interpolation.__name__}: tvec is not \
                         sorted or have values < 0.0 !\nCheck sheet {sheet} \
                         in {f_path}\n"
            )
    return yy


def read_interp_file(file_path, comp, *INTIAL, **options):

    """
    ##############################################################################
    #              Read_interp_file(file_path, comp)
    ##############################################################################
    #
    # Function that reads data to be interpolated from .xlsx files.
    #
    ##############################################################################
    # VARIABLE    I/O    TYPE                  DESCRIPTION                  UNIT
    # --------------------------------------------------------------------------
    # file_path   I      string                path of the file to be read  -
    # comp        I      object                python object of
    #                                          class SolidComponent    -
    # tt          O      np array float        time values vector to be
    #                                          used in interpolation        s
    # xx          O      np array float        space values vector to be
    #                                          used in interpolation        m
    # Data_matrix O      np 2 dimensional      matrix of data to be
    #                    array float           used in interpolation        many
    # flag        O      scalar integer        flag that specifies
    #                                          data units                   -
    ##############################################################################
    #
    # Invoked functions/methods: none
    #
    ##############################################################################
    # N.B. no self parameters in input.
    # N.B. In case of bfield.xlsx and alphaB.xlsx if flag == 1 data unit are T
    # and T/m respectively, if flag == 2 data are per current value (i.e T/A
    # and T/m/A)
    ##############################################################################
    #
    # Author D. Placido Polito 06/2020
    #
    ##############################################################################
    """

    # Function Read_interp_file starts here. (cdp)

    ff = load_workbook(file_path, data_only=True)
    sheet = ff[comp.ID]
    flag = sheet.cell(1, 1).value  # read the element in cell (1,1) of .xlsx file
    # get sheet boundary to read data properly (cdp, 07/2020)
    r_first = sheet.min_row + 1  # first row to be read (cdp, 07/2020)
    r_last = sheet.max_row  # last row to be read (cdp, 07/2020)
    c_first = sheet.min_column  # first column to be read (cdp, 07/2020)
    c_last = sheet.max_column  # last column to be read (cdp, 07/2020)
    r_start = r_first + 1

    if options.get("action") == "time_&_space":
        # interpolation in both time and space (cdp, 07/2020)
        c_start = c_first + 1
        # Data matrix initialization
        Data_matrix = np.zeros((r_last - r_first, c_last - c_first))
        # for loop to read tt values form file
        for cc in sheet.iter_rows(
            min_row=r_first,
            max_row=r_first,
            min_col=c_start,
            max_col=c_last,
            values_only=True,
        ):
            # get time @ which data are given in a single shot (cdp, 07/2020)
            tt = np.array(cc)
        # for loop to read xx values form file
        for rr in sheet.iter_cols(
            min_row=r_start,
            max_row=r_last,
            min_col=c_first,
            max_col=c_first,
            values_only=True,
        ):
            # get spatial coordinates @ which data are given in a single shot \
            # (cdp, 07/2020)
            xx = np.array(rr)
        # Data_matrix construction
        for rr in range(r_start, r_last + 1):
            for cc in sheet.iter_rows(
                min_row=rr,
                max_row=rr,
                min_col=c_start,
                max_col=c_last,
                values_only=True,
            ):
                Data_matrix[rr - r_start, :] = np.array(cc)
        return [flag, tt, Data_matrix, xx, sheet]
    elif options.get("action") == "time":
        # interpolation in time, only for flow variables at channel inlet \
        # (cdp, 07/2020)
        c_start = c_first + 2
        # for loop to read tt values form file
        for cc in sheet.iter_rows(
            min_row=r_first,
            max_row=r_first,
            min_col=c_start,
            max_col=c_last,
            values_only=True,
        ):
            # get time @ which data are given in a single shot (cdp, 07/2020)
            tt = np.array(cc)
        # read input data names
        for rr in sheet.iter_cols(
            min_row=r_start,
            max_row=r_last,
            min_col=c_first,
            max_col=c_first,
            values_only=True,
        ):
            val_name = np.array(rr, dtype=str)
        if INTIAL[0] == -1:
            # get TEMINL, TEMOUT, PREINL and PREOUT rows (cdp, 07/2020)
            req_variable = np.array(["TEMINL", "TEMOUT", "PREINL", "PREOUT"], dtype=str)
            row_ind = np.zeros(req_variable.shape, dtype=int)
            for ii in range(len(req_variable)):
                row_ind[ii] = int(np.nonzero(val_name == req_variable[ii])[0] + r_start)
        elif INTIAL[0] == -2:
            # get TEMINL, TEMOUT, PREINL and MDTIN rows (cdp, 07/2020)
            req_variable = np.array(["TEMINL", "TEMOUT", "PREINL", "MDTIN"], dtype=str)
            row_ind = np.zeros(req_variable.shape, dtype=int)
            for ii in range(len(req_variable)):
                row_ind[ii] = int(np.nonzero(val_name == req_variable[ii])[0] + r_start)
        elif INTIAL[0] == -5:
            # get TEMINL, TEMOUT, PREOUT ans MDTIN rows (cdp, 07/2020)
            req_variable = np.array(["TEMINL", "TEMOUT", "PREOUT", "MDTIN"], dtype=str)
            row_ind = np.zeros(req_variable.shape, dtype=int)
            for ii in range(len(req_variable)):
                row_ind[ii] = int(np.nonzero(val_name == req_variable[ii])[0] + r_start)

        # create Data_matrix by row: each row correspond to a different input \
        # values according to the above if else series. (cdp, 07/2020)

        # Data_matrix initialization (cdp, 07/2020)
        Data_matrix = np.zeros((row_ind.shape[0], tt.shape[0]))
        for ii in range(len(row_ind)):
            for values in sheet.iter_rows(
                min_row=row_ind[ii],
                max_row=row_ind[ii],
                min_col=c_start,
                max_col=c_last,
                values_only=True,
            ):
                Data_matrix[ii, :] = np.array(values)
        return [flag, tt, Data_matrix, sheet]


# Function pt_binsrc.py starts here
def pressure_temperature_binary_search(P, T, PHECUR, IPHECUR, dict_He_tables, NPRESS):

    """
    ######################################################################
    #
    # Binary search in pressure and temperature.
    # Returns PERFECT1/2 if the T is greater then the last T.
    #
    # Author : Bonifetto R., Energy Department, Politecnico di Torino.
    # Version: 1.0 (July 2012).
    # Translation form MatLab to Python: D.Placido Polito 25/03/2020
    # Correction after test: D.Placido PoliTo 10/04/2020
    ######################################################################
    """
    if P.shape != T.shape:
        raise ValueError(
            f"ERROR: pressure and temperature array must have the same shape:\nP_shape = {P.shape}\nT_shape = {T.shape}\n"
        )
    # End if P.shape

    # Variables initialization
    IP = -1 * np.ones(P.shape, dtype=int)
    IT1 = -1 * np.ones(P.shape, dtype=int)
    F1 = np.zeros(P.shape)
    PERFECT1 = -1 * np.ones(P.shape, dtype=int)
    TT1 = -1 * np.ones(P.shape)
    F2 = np.zeros(P.shape)
    IT2 = -1 * np.ones(P.shape, dtype=int)
    PERFECT2 = -1 * np.ones(P.shape, dtype=int)
    TT2 = -1 * np.ones(P.shape)

    PP = P

    # Placido Daniele comment: it seems that HEPRESS and HETEMP are arrays while \
    # IP, IT1 and IT2 are some index to access these vectors. Since this is a \
    # translation from MatLab (starts counting from 1) to Python (starts \
    # counting from 0), I decrease by 1 all the index to be consistent.
    # NPRESS is number of data in Pressure, typically fixed at 290. Since \
    # counting starts from 0, to have NPRESS (290) number consistent it becomes \
    # NPRESS-1
    # This function will return IP = IP - 1, IT1 = IT1 - 1 and IT2 = IT2 - 1 \
    # instead of IP, IT1 and IT2 MatLab values. (26/03/2020)

    NPRESS = NPRESS - 1
    # Check range of Pressure
    ind = np.nonzero(
        (PP < dict_He_tables["HEPRESS"][0]) | (PP > dict_He_tables["HEPRESS"][-1])
    )[0]

    PP[ind] = np.maximum(dict_He_tables["HEPRESS"][0], PP[ind])
    PP[ind] = np.minimum(dict_He_tables["HEPRESS"][-1], PP[ind])

    for ii in range(len(P)):

        # Same pressure ?
        if (PP[ii] == PHECUR[ii]) and (
            IPHECUR[ii] != 0
        ):  # this condition is slightly different from the one in HeRef.py (Placido Daniele, 28/03/2020)
            IP[ii] = IPHECUR[ii] - 1
        else:
            # Binary search on Pressure
            IP[ii] = bisect.bisect_left(dict_He_tables["HEPRESS"], PP[ii]) - 1
            PHECUR[ii] = PP[ii]
            IPHECUR[ii] = IP[ii]

        # Curve at dict_He_tables["HEPRESS"][IP], get number of data points and \
        # positions
        IND1S = 0  # IND1S = 1 (old MatLab) put equal 0 since python starts \
        # counting from 0 and not from 1 as MatLab does (Placido Daniele 26/03/2020)
        IND1E = dict_He_tables["NTEMP"][IP[ii]] - 1

        # Curve at dict_He_tables["HEPRESS"][IP + 1], get number of data points \
        # and positions
        IND2S = 0  # IND2S = 1 (old MatLab) put equal 0 since python starts \
        # counting from 0 and not from 1 as MatLab does (Placido Daniele 26/03/2020)
        IND2E = dict_He_tables["NTEMP"][IP[ii] + 1] - 1

        # Perfect gas condition: no interpolation
        PERFECT1[ii] = T[ii] >= dict_He_tables["HETEMP"][IP[ii], IND1E]
        PERFECT2[ii] = T[ii] >= dict_He_tables["HETEMP"][IP[ii] + 1, IND2E]

        # Binary search in temperature
        # Not a perfect gas (cdp, 08/2020)
        if PERFECT1[ii] == 0:
            if T[ii] < dict_He_tables["HETEMP"][IP[ii], IND1S]:
                # This is consistent with what is done in 4C \
                # (see PT_BINSRC HeRefProp.f90) (cdp, 10/2020)
                warnings.warn(
                    "Saturation 1: IT1 forced to 0 to avoid code crushing and TT1 set equal to the saturation temperature at given pressure.\n"
                )
                IT1[ii] = 0
                TT1[ii] = dict_He_tables["HETEMP"][IP[ii], IND1S]
            else:
                TT1[ii] = max(T[ii], dict_He_tables["HETEMP"][IP[ii], IND1S])
                IT1[ii] = (
                    bisect.bisect_left(
                        dict_He_tables["HETEMP"][IP[ii], IND1S:IND1E:1], TT1[ii]
                    )
                    - 1
                )
            # End if T[ii] (cdp, 10/2020)
        else:
            TT1[ii] = T[ii]
            IT1[ii] = IND1E
        # End if PERFECT1[ii] (cdp, 10/2020)

        # Not a perfect gas (cdp, 08/2020)
        if PERFECT2[ii] == 0:
            if T[ii] < dict_He_tables["HETEMP"][IP[ii] + 1, IND2S]:
                # This is consistent with what is done in 4C \
                # (see PT_BINSRC HeRefProp.f90) (cdp, 10/2020)
                warnings.warn(
                    "Saturation 2: IT2 forced to 0 to avoid code crushing TT1 set equal to the saturation temperature at a pressure value in table that is in the row below the one of the given pressure.\n"
                )
                IT2[ii] = 0
                TT2[ii] = dict_He_tables["HETEMP"][IP[ii] + 1, IND2S]
            else:
                TT2[ii] = max(T[ii], dict_He_tables["HETEMP"][IP[ii] + 1, IND2S])
                IT2[ii] = (
                    bisect.bisect_left(
                        dict_He_tables["HETEMP"][IP[ii] + 1, IND2S:IND2E:1], TT2[ii]
                    )
                    - 1
                )
            # End if T[ii] (cdp, 10/2020)
        else:
            TT2[ii] = T[ii]
            IT2[ii] = IND2E
        # End if PERFECT2[ii] (cdp, 10/2020)
    # End for loop

    # array smart index (cdp, 08/2020)
    ind = np.nonzero(PERFECT1 == 0)[0]
    # Used only if He is not a perfect gas (cdp, 08/2020)
    F1[ind] = (TT1[ind] - dict_He_tables["HETEMP"][IP[ind], IT1[ind]]) / (
        dict_He_tables["HETEMP"][IP[ind], IT1[ind] + 1]
        - dict_He_tables["HETEMP"][IP[ind], IT1[ind]]
    )

    # Overwriting ind (cdp, 06/2020)
    ind = np.nonzero(PERFECT2 == 0)[0]
    # Used only if He is not a perfect gas (cdp, 08/2020)
    F2[ind] = (TT2[ind] - dict_He_tables["HETEMP"][IP[ind] + 1, IT2[ind]]) / (
        dict_He_tables["HETEMP"][IP[ind] + 1, IT2[ind] + 1]
        - dict_He_tables["HETEMP"][IP[ind] + 1, IT2[ind]]
    )

    return [PP, PHECUR, IPHECUR, IP, IT1, F1, PERFECT1, IT2, F2, PERFECT2]


# End function PT_BINSRC.

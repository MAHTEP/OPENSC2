import numpy as np


class CylindricalHelix:
    """Class with methods and attributes to evaluate some cylindrical helix geometrical information."""

    def __init__(self, x0: float, y0: float, heigth: float, costheta: float):
        """Method that make instances of class Helix. Allow to define an helix with some of the main geometrical information such as:
        * radius;
        * pitch angle;
        * reduced pitch;
        * number of windings
        * length of a single winding
        * total helix length

        Args:
            x0 (float): abscissa of the barycenter of the generic  StrandMixedComponent or StrandStabilizerComponent or StrandSuperconductorComponent object.
            y0 (float): ordinate of the barycenter of the generic  StrandMixedComponent or StrandStabilizerComponent or StrandSuperconductorComponent object.
            heigth (float): heigth of the cylinder around which the helix wounds.
            costheta (float): cosine of the pitch angle of the helix.
        """
        # Height of the helix.
        self.height = heigth
        # Evaluate cylinder radius
        self.__eval_radius(x0, y0)
        # Evaluate angle in the xy plane.
        self.__eval_alpha(x0, y0)
        # Evaluate pitch angle.
        self.__eval_theta(costheta)
        # Evaluate reduced_pitch.
        self.__eval_reduced_pitch()
        # Evaluate windings number.
        self.__eval_winding_number()
        # Evaluate the length of a single winding.
        self.__eval_winding_length()
        # Evaluate the total lenght of the helix.
        self.__eval_helix_length()

    def __eval_radius(self, x0: float, y0: float):
        """Private method that evaluates the radius of the cylinder around which the helix wounds starting from x0 and y0 coordinates (the coorinates of the barycenter of the generic StrandMixedComponent or StrandStabilizerComponent or StrandSuperconductorComponent object). Value is assigned to attribute self.radius.

        Args:
            x0 (float): abscissa of the barycenter of the generic  StrandMixedComponent or StrandStabilizerComponent or StrandSuperconductorComponent object.
            y0 (float): ordinate of the barycenter of the generic  StrandMixedComponent or StrandStabilizerComponent or StrandSuperconductorComponent object.
        """
        self.radius = np.sqrt(x0 ** 2 + y0 ** 2)

    def __eval_alpha(self, x0: float, y0: float):
        """Private method that evaluates the angle in the xy plane; used in helix parametrization. Value is assigned to attribute self.alpha.

        Args:
            x0 (float): abscissa of the barycenter of the generic  StrandMixedComponent or StrandStabilizerComponent or StrandSuperconductorComponent object.
            y0 (float): ordinate of the barycenter of the generic  StrandMixedComponent or StrandStabilizerComponent or StrandSuperconductorComponent object.
        """
        self.alpha = np.arctan2(y0, x0)

    def __eval_theta(self, costheta: float):
        """Private method that evaluates the pitch angle of the helix. Value is assigned to attribute self.theta.

        Args:
            costheta (float): cosin of the pitch angle.
        """
        self.theta = np.arccos(costheta)

    def __eval_reduced_pitch(self):
        """Private method that evaluates the reduced pitch of the helix. Value is assigned to attribute self.reduced_pitch."""
        # From geometrical definition of tan(theta).
        self.reduced_pitch = self.radius / np.tan(self.theta)

    def __eval_winding_number(self):
        """Private method that evaluates the number of windings of the helix around its heligth. Value is assigned to attribute self.winding_number."""
        self.winding_number = self.height / (2 * np.pi * self.reduced_pitch)

    def __eval_winding_length(self):
        """Private method that evaluates the length of a single winding of the helix. Value is assigned to attribute self.winding_length."""
        self.winding_length = (
            2 * np.pi * np.sqrt(self.radius ** 2 + self.reduced_pitch ** 2)
        )

    def __eval_helix_length(self):
        """Private method that evaluates the total length of the helix. Value is assigned to attribute self.length."""
        self.length = self.winding_number * self.winding_length

    def helix_parametrization(self, tau: np.ndarray) -> tuple:
        """Method that evaluates cylindrical helix coordinates according to canonical parametrization.

        Args:
            tau (np.ndarray): angular discretization.

        Returns:
            tuple: x, y and z coordinate of the cylindrical helix.
        """
        xx = self.radius * np.cos(self.alpha + tau)
        yy = self.radius * np.sin(self.alpha + tau)
        zz = self.reduced_pitch * tau
        return xx, yy, zz
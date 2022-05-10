class ComponentCollection:
    """Class that allows to organize the components that build the conductor in collection(s).
    A collection is characterized by:
        1) a list of objects (not necessary of the same kind);\n
        2) the number of objects that make up the collection;\n
        3) the collection name (optional parameter).
    """

    def __init__(self, name=None):
        """Make an istance of class ComponentCollection.

        Args:
            name (str, optional): name of the collection. Defaults to None.
        """
        self.collection = list()
        self.number = len(self.collection)
        self.name = name
class ComponentsInventory:
    """Class that allows to organize the components that build the conductor """

    def __init__(self, name=None):
        """Make an istance of class ComponentsInventory.

        Args:
            name (str, optional): name of the collection. Defaults to None.
        """
        self.collection = list()
        self.number = len(self.collection)
        self.name = name
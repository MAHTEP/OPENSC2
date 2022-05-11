import numpy as np

class CylindricalHelix:
    """Class with methods and attributes to evaluate some cylindrical helix geometrical information."""

    def __init__(self, radius: float, heigth: float, costheta: float):
        """Method that make instances of class Helix. Allow to define an helix with some of the main geometrical information such as:
        * reduced pitch;
        * number of windings
        * length of a single winding
        * total helix length

        Args:
            radius (float): radius of the cylinder around which the helix wounds.
            heigth (float): heigth of the cylinder around which the helix wounds.
            costheta (float): cosine of the pitch angle of the helix.
        """
        # Radius of the helix.
        self.radius = radius
        # Height of the helix.
        self.height = heigth
        # Evaluate theta.
        self.__eval_theta(costheta)
        # Evaluate reduced_pitch.
        self.__eval_reduced_pitch()
        # Evaluate windings number.
        self.__eval_winding_number()
        # Evaluate the length of a single winding.
        self.__eval_winding_length()
        # Evaluate the total lenght of the helix.
        self.__eval_helix_length()


    def __eval_theta(self, costheta: float):
        """Private method that evaluate the pitch angle of the helix. Value is assigned to attribute self.theta.

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
        """Private method that evaluate the length of a single winding of the helix. Value is assigned to attribute self.winding_length."""
        self.winding_length = (
            2 * np.pi * np.sqrt(self.radius ** 2 + self.reduced_pitch ** 2)
        )

    def __eval_helix_length(self):
        """Private method that evaluate the total length of the helix. Value is assigned to attribute self.length."""
        self.length = self.windings_number * self.winding_length

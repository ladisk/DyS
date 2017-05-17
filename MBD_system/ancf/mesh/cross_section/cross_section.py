"""

created by: lskrinjar
date of creation: 03/08/2016
time of creation: 16:24
"""
import numpy as np

class CrossSection(object):
    """
    Mesh object constructor
    """

    def __init__(self, type=None, parent=None):
        """

        :param type:
        :param parent:
        """
        #   parent
        self._parent = parent
        #   type of cross section
        self.type = type
        #   cross section supported types
        #   rectangle
        #   rectangular tube
        #   circular
        #   circular tube

        #   area
        self.area = None
        #   inertia of area
        self.Iyy = None
        self.Izz = None

        #   rectangle
        #   width
        self.B = None
        #   height
        self.H = None

        #   circular
        self.R = None

        self.__info_cross_section_not_defined = "Cross section geometry parameters are not defined! Evaluation nod executed!"

    def evaluate(self):
        """

        :return:
        """
        if self.type == "rectangle":
            if self.B is not None and self.H is not None:
                self._evaluate_rectangle()
            else:
                print self.__info_cross_section_not_defined

        elif self.type == "circle":
            if self.R is not None:
                self._evaluate_circle()
            else:
                print self.__info_cross_section_not_defined

        else:
            raise ValueError, "Cross section type not correct!"

    def _evaluate_rectangle(self):
        """

        :return:
        """
        self.area = CrossSection.area_rectangle(self.B, self.H)

        self.Izz = CrossSection.I_rectangle(self.B, self.H)

        self.Iyy = CrossSection.I_rectangle(self.H, self.B)

    def _evaluate_circle(self):
        """

        :return:
        """
        self.area = CrossSection.area_circular(self.R)

        self.Izz = self.Iyy = CrossSection.I_circular(self.R)

    @staticmethod
    def area_rectangle(b, h):
        return b * h

    @staticmethod
    def I_rectangle(b, h):
        return (b * h**3) / 12.

    @staticmethod
    def area_circular(R):
        return np.pi * R**2

    @staticmethod
    def I_circular(R):
        return (np.pi / 4.) * (R**4)


if __name__ == "__main__":
    section = CrossSection(type="circle")
    section.R = 0.057406263
    section.evaluate()
    print "A =", section.area
    print "Izz =", section.Izz


"""

created by: lskrinjar
date of creation: 16/03/2016
time of creation: 16:34
"""
import itertools
import numpy as np


from body import Body
from grid.grid import Grid


class GroundBody(Body):
    """
    classdocs
    """

    def __init__(self, parent=None):
        """

        :param parent:
        :return:
        """
        super(GroundBody, self).__init__(name="ground", parent=parent)

        #   type of body
        self.body_type = "ground"

        #   body id
        self.body_id = -1

        #   name
        self._name = "Ground"

        #   visualization properties
        self.color = np.array([0.8, 0.8, 0.8])

        #   grid object
        self.grid = Grid(parent=self)

    def set_vtk_data(self):
        """
        Create grid to display
        :return:
        """
        self.grid.set_vtk_data()

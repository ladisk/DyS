"""

created by: lskrinjar
date of creation: 04/03/2016
time of creation: 17:17
"""
import numpy as np
from OpenGL.GL import *


from geometry import Geometry


class Geometry2D(Geometry):
    """
    classdocs
    """

    def __init__(self, filename, parent=None):
        """

        :param filename:
        :param parent:
        :return:
        """
        super(Geometry2D, self).__init__(filename, parent=parent)



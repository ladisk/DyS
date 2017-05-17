"""
Created on 10. jan. 2015

@author: lskrinjar
"""
import numpy as np

class DataContainerAABB(object):
    """
    classdocs
    """

    def __init__(self, _id=None, _nodes=[], _normals=[], _tangents=[], _angles=[], x_min=None, x_max=None, y_min=None, y_max=None, _type=None):
        """
        Constructor of object that saves data of sub AABB to be used in next time step in parent in recursive built of AABB tree
        Args:
            _id         - id number of an object
            _nodes      - nodes (a matrix - ndarray) of nodes
            _normals    - normals (a matrix - ndarray) of normals
            x_min       - min x coordinate of nodes (float)
            x_max       - max x coordinate of nodes (float)
            y_min       - min y coordinate of nodes (float)
            y_max       - max y coordinate of nodes (float)
        """
        self._id = _id
        #   type is direction of division and min or max subsection
        self._type = _type

        self._nodes = _nodes
        self._normals = _normals
        self._tangents = _tangents
        self._angles_LCS = _angles

        #   dimensions
        self.__dims(x_min, x_max, y_min, y_max)

    def __dims(self, x_min, x_max, y_min, y_max):
        """
        
        """
        #    min and max dimensions
        self._min = np.array([x_min, y_min])
        
        self._max = np.array([x_max, y_max])
                
        self._constructed = True
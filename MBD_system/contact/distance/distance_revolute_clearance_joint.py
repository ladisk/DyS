__author__ = 'lskrinjar'

import numpy as np

class DistanceRevoluteClearanceJoint(object):
    """

    """
    def __init__(self, u_iP, u_jP, parent=None):
        """
        Constructor of a distance object
        :param u_iP:    center of body i hole in GCS coordinates
        _param u_jP:    center of body j pin in GCS coordinates
        :return:
        """
        #   parent
        self._parent = parent

        self._distance_vector = u_jP - u_iP

        self._distance = np.linalg.norm(self._distance_vector, ord=2)
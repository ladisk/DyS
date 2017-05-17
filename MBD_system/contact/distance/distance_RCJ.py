__author__ = 'lskrinjar'

import numpy as np
from MBD_system.n2t import n2t


class DistanceRCJ(object):
    """
    classdocs
    """

    def __init__(self, r_iP, r_jP, parent=None):
        """
        Constructor of a distance object
        :param u_iP:    center of body i hole in GCS coordinates
        _param u_jP:    center of body j pin in GCS coordinates
        :return:
        """
        #   parent
        self._parent = parent

        #   body ids
        #   i - free point of pin body id
        #   j - line body of slot body id
        self.body_id_i = self._parent.body_id_i
        self.body_id_j = self._parent.body_id_j
        self.body_id_list = [self.body_id_i, self.body_id_j]

        #   points in GCS
        self.r_jP = r_jP
        self.r_iP = r_iP

        self._distance_vector = self.r_jP - self.r_iP

        self._distance = np.linalg.norm(self._distance_vector, ord=2)

        #    normal in GCS of cylindrical surface
        self.normal = self._distance_vector / self._distance

        #    tangent in GCS
        self.tangent = n2t(self.normal)



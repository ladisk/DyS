__author__ = 'lskrinjar'

import numpy as np
from matplotlib import pyplot as plt

from MBD_system.A import A_matrix


class DistanceLineNode(object):#Distance
    """
    classdocs
    """
    def __init__(self, u_iP, u_jP, normal=None, u_iR=None, tangent=None, parent=None):
        # super(DistancePlaneNode, self).__init__(self._name, parent)
        """
        :param u_iP:        free a point in GCS on body i
        :param u_jP:        line point in GCS on body j
        :param normal:      normal in GCS
        :param u_iR:        line point in GCS on body j
        :param tangent:     tangent in GCS
        """
        #   parent
        self._parent = parent

        #   coordinates in GCS
        #   body i - edge
        #   body j - node
        self.u_iP = u_iP
        self.u_jP = u_jP
        self.u_iR = u_iR

        #   edge vector
        if u_iR is not None:
            self.edge = u_iR - u_iP
        else:
            self.edge = None

        #   distance
        self._distance_vector = u_jP - u_iP

        #   tangent
        if tangent is None and normal is None:
            self.tangent = self.edge / np.linalg.norm(self.edge, ord=2)
        elif tangent is None and normal is not None:
            self.tangent = np.dot(A_matrix(np.pi/2), normal)
        else:
            self.tangent = tangent

        #   normal
        if normal is None:
            self.normal = np.array([self.tangent[1], -self.tangent[0]])
        else:
            self.normal = normal

        #   distance projection
        _distance_vector_projection = self.normal * np.dot(self._distance_vector, self.normal)

        #   distance
        self._distance = np.linalg.norm(_distance_vector_projection, ord=2)
        # print "self._distance =", self._distance

        #    distance sign
        # print "np.cross(self.tangent, _distance_vector) =", np.cross(self.tangent, _distance_vector)
        if np.sign(np.cross(self.tangent, self._distance_vector)) < 0:
            self._inside = False
            self._distance_sign = self._distance
        else:
            self._inside = True
            self._distance_sign = -self._distance

    def _evaluate_CP_GCS(self):
        """

        :return:
        """
        self.CP = self.u_iP + np.dot(self._distance_vector, self.tangent) * self.tangent

    def contact_point_on_line_GCS(self):
        """
        Function calculates contact point - CP
        :return:
        """
        self._evaluate_CP_GCS()
        return self.CP


if __name__ == "__main__":
    #   free node
    u_iP = np.array([1.5, -1.1346])
    #   line node
    u_jP = np.array([0, 0])
    normal = np.array([0, 1])
    d = DistanceLineNode(u_iP, u_jP, normal)
    print d._distance

    fig = plt.figure(num=1, figsize=(6, 4), dpi=100, facecolor='w', edgecolor='k')
    ax = plt.subplot(111, aspect="auto")
    #   normal
    plt.quiver([normal[0], normal[0]], [0, normal[1]], linestyle="-", color=(0, 0, 0))
    #   line node
    plt.plot(u_jP[0], u_jP[1], marker="x", color=(0, 0, 0))
    #   free node
    plt.plot(u_iP[0], u_iP[1], marker="x", color=(1, 0, 0), markersize=5)
    #   distance
    plt.plot([u_iP[0], u_iP[0]], [0, d._distance], linestyle="-", color=(0, 1, 0), markersize=5)

    plt.show()
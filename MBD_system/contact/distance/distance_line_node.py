__author__ = 'lskrinjar'

from pprint import pprint
import numpy as np
from matplotlib import pyplot as plt

from MBD_system.A import A_matrix
from MBD_system.contact.distance.distance import Distance

class DistanceLineNode(object):#Distance
    """
    classdocs
    """
    def __init__(self, u_iP, u_jP, normal, u_iR=None, tangent=None, parent=None):
        # super(DistancePlaneNode, self).__init__(self._name, parent)
        """
        :param u_iP:        free a point in GCS
        :param u_jP:        line point in GCS
        :param normal:      normal in GCS
        :param u_iR:        a point in GCS
        :param tangent:     tangent in GCS
        """
        #   parent
        self._parent = parent

        #   distance
        _distance_vector = u_iP - u_jP

        #   tangent
        self.tangent = np.dot(A_matrix(np.pi/2), normal)

        _distance_vector_projection = normal * np.dot(_distance_vector, normal)

        self._distance = np.linalg.norm(_distance_vector_projection, ord=2)
        # print "self._distance =", self._distance

        #    distance sign
        # print "np.cross(self.tangent, _distance_vector) =", np.cross(self.tangent, _distance_vector)
        if np.sign(np.cross(self.tangent, _distance_vector)) < 0:
            self._inside = False
            self._distance_sign = self._distance
        else:
            self._inside = True
            self._distance_sign = -self._distance

        # print "self._inside =", self._inside

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
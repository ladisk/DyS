'''
Created on 19. mar. 2015

@author: lskrinjar
'''
from pprint import pprint
import numpy as np
from matplotlib import pyplot as plt

class Distance(object):
    """
    classdocs
    """
    def __init__(self, n0, n1, n2=None, normal=None, node_body=None, edge_body=None, parent = None):
        """
        Constructor
        """
        self._parent = parent

        self._distance = np.inf

        self.node = n0
        self.edge_node_1 = n1
        self.edge_node_2 = n2

        self.normal = normal
        self.tangent = None

        self.node_body = node_body
        self.edge_body = edge_body

        #   max penetration depth attribute assigned from body attribute
        if edge_body is not None:
            self._max_penetration_depth = edge_body._max_penetration_depth
        else:
            self._max_penetration_depth = 1E-3

        #   evaluate distance
        if n2 is None:
            self._distance = self.__vector_length(n0, n1)
        else:
            self.__evaluate_distance_2D(n0, n1, n2)

    def __evaluate_distance_2D(self, n0, n1, n2):
        """
        Args:
            n0 - free node
            n1 - start node of edge
            n2 - end node of edge 
        """
        #   vector from start node of edge to free node
        n1n0 = n1 - n0
        n2n0 = n2 - n0
        #   edge vector from 2 to 1
        n2n1 = n2 - n1

        #   calculate angles
        #   angle
        fi_012 = self.angle(n1n0, -n2n1)
        fi_021 = self.angle(n2n0, n2n1)

        #   tangent
        self.__tangent = n2n1
        self.tangent = n2n1_e = n2n1 / np.linalg.norm(n2n1)

        _distance_vector = (np.dot(n1n0, n2n1_e) * n2n1_e) - n1n0

        self._distance = np.linalg.norm(_distance_vector, ord=2)

        normal_plus_distance = self.normal + _distance_vector

        #   check if free node is inside (TF status)
        if (np.linalg.norm(normal_plus_distance) < np.linalg.norm(self.normal)) and (self._distance <= self._max_penetration_depth) and (fi_012 <= np.pi/2) and (fi_021 <= np.pi/2):
            self._inside = True
            self._distance_sign = -self._distance
        else:
            self._inside = False
            self._distance_sign = +self._distance

    def __vector_length(self, T1, T2):
        """
        Vector length is calculated as vector norm
        """
        #   vector formed between points T1, T2
        _distance_vector = T2 - T1

        #   vector length
        _d = np.linalg.norm(_distance_vector, ord=2)

        #   calculate normal for revolute clearance joint
        self.normal = _distance_vector/_d

        self.tanget = None
        return _d

    def angle(self, vec1, vec2):
        return np.arccos(np.dot(vec1, vec2) / (np.linalg.norm(vec1, ord=2) * np.linalg.norm(vec2, ord=2)))

    def get_normal_2D(self):
        return self.normal[0:2]

    def get_tangent_2D(self):
        return self.tangent[0:2]

    def get_edge_2D(self):
        return np.array([self.edge_node_1[0:2], self.edge_node_2[0:2]])

    def get_node_2D(self):
        return np.array(self.node[0:2])

    def distance_2D_in_direction_of_normal(self):
        """
        
        """
        self._distance_2D = 0

    def set_data_LCS(self, n0_LCS, n1_LCS, n2_LCS):
        """
        Function calculates contact data in LCS of each body in contact pair (node, edge)

        :return:
        """


if __name__ == "__main__":
    n0 = np.array([.08, .10, 0])
    n1 = np.array([0.1, 0.1, 0])
    n2 = np.array([.1, 0, 0])
    normal = np.array([1, 0, 0])
    # n0 = np.array([.08, .1, 0])
    # n1 = np.array([0.1, 0.05, 0])
    # n2 = np.array([.1, 0.1, 0])


    d = Distance(n0, n1, n2, normal)
    pprint(vars(d))


    fig = plt.figure(num=0, figsize=(6, 4), dpi=80, facecolor='w', edgecolor='k')
    ax = fig.add_subplot(111)
    plt.grid()
    plt.xlim([0.06, +0.12])
    plt.ylim([-0.02, +0.12])

    #   free node
    ax.plot(n0[0], n0[1], color="black", marker="o", linestyle = "-")
    ax.text(n0[0]*1.01, n0[1]*1.01, "N0 =["+str(n0[0])+","+str(n0[1])+"]", color='black')

    #   edge
    ax.arrow(n1[0], n1[1], n2[0]-n1[0], n2[1]-n1[1], width = 1E-4, head_width=4E-4, head_length=0.01, fc='k', ec='k')
    ax.plot([n1[0], n2[0]], [n1[1], n2[1]], color="black", marker="s", linestyle = "-")
    ax.text(n1[0]*1.01, n1[1]*1.01, "N1 =["+str(n1[0])+","+str(n1[1])+"]", color='black')
    ax.text(n2[0]*1.01, n2[1]*1.01, "N2 =["+str(n2[0])+","+str(n2[1])+"]", color='black')
    #   normal
    ax.arrow(n1[0], n1[1], normal[0]*1E-2, normal[1]*1E-2, width = 1E-4, head_width=4E-4, head_length=1E-4, fc='grey', ec='grey')

    #   plot n0n1
    ax.plot([n0[0], n1[0]], [n0[1], n1[1]], color="grey", linestyle = "--")

    plt.show()
"""
Created on 19. mar. 2015

@author: lskrinjar
"""
from pprint import pprint

import numpy as np
from matplotlib import pyplot as plt

from MBD_system.Ai_ui_P import Ai_ui_P_vector
from MBD_system.q2R_i import q2R_i
from MBD_system.q2theta_i import q2theta_i
from MBD_system.n2t import n2t


class Distance(object):
    """
    classdocs
    """
    def __init__(self, r_iP, r_jP, r_jR=None, normal=None, q=None, node_body=None, edge_body=None, parent=None):
        """
        Constructor of a distance object that calculates distance (value, scalar) of free node to edge
        :param n0:          free node of a body in body LCS
        :param r_jP:          start node of an edge of a body in body LCS
        :param r_jR:          end node of an edge of a body in body LCS
        :param normal:      normal of an edge in body LCS
        :param q:           vector of generalized coordinates at time t
        :param node_body:   pointer to node body object
        :param edge_body:
        """
        #   parent
        self._parent = parent

        #   default value for attribute
        self._distance = np.inf

        #   initialize data containers
        self._data_container()

        #    node on body i - free node
        self.r_iP = r_iP

        #   node on body j
        self.r_jP = r_jP

        #    normal
        self.normal = normal

        #   tangent
        if self.normal is not None:
            self.tangent = n2t(self.normal)
        else:
            self.tangent = None

        if node_body is not None and edge_body is not None:
            #    pointers to edge and node body
            self.node_body = node_body
            self.node_body_id = node_body.body_id
            self.edge_body = edge_body
            self.edge_body_id = edge_body.body_id

            #   edge vector from 1 to 2
            self.r_jRr_jP = r_jR - r_jP
            # print "self.edge =", self.r_jRr_jP
            # angle = np.arccos(np.dot(r_jP, r_jR))
            # print "angle =", angle
            
            #    edge nodes
            self.r_jP = r_jP
            self.r_jR = r_jR
            #    edge vector
            self.edge = r_jR - r_jP

            #    evaluate direction of edge vector it has to be in direction of normal
            self._check_direction(r_jP, r_jR)
    
            #   unit edge vector
            self.r_jRr_jP_e = self.edge / np.linalg.norm(self.edge)
                
        #   max penetration depth attribute assigned from body attribute
        if edge_body is not None:
            self.max_penetration_depth = edge_body.max_penetration_depth
        else:
            self.max_penetration_depth = 1E-4
        
        self.min_penetration_depth = 10*self.max_penetration_depth

        #   evaluate distance
        #   evaluate distance: node to node
        if r_jR is None:
            self._distance, self._inside = self._vector_length(self.r_iP, self.r_jP)
        #   evaluate distance: node to edge
        else:
            self._evaluate_distance_2D(q)
        
        #    calculate tangent based on normal
        self.tangent = Ai_ui_P_vector(self.normal, np.pi/2)
            
    def _check_direction(self, r_jP, r_jR):
        """
        Function checks direction of edge vector.
        If direction of edge vector is equal as direction of tangent (calculated from normal) it is OK, else the values
        of attributes r_jP and r_jR are replaced
        r_jP -> r_jR
        r_jR -> r_jP
        """
        print "self.tangent =", self.tangent
        print "self.edge =", self.edge
        if (np.sign(self.tangent) == np.sign(self.edge)).all():
            pass
        else:
            self.r_jP = r_jR
            self.r_jR = r_jP
            #    edge vector
            self.edge = r_jP - r_jR

    def _data_container(self):
        """
        Function initializes container to store contact data at each time step, such as
        distance value, time, inside status
        :return:
        """
        self._t_solution_container = [0]
        self._distance_solution_container = [0]

    def _evaluate_distance_2D(self, q=None):
        """
        Args:
            n0 - free node
            r_jP - start node of edge
            r_jR - end node of edge 
        """
        #   coordinates in GCS
        # print "self.n0 =", self.n0
        # print "q2R_i(q, self.node_body_id) =", q2R_i(q, self.node_body_id)
        # print "q2theta_i(q, self.node_body_id) =", q2theta_i(q, self.node_body_id)
        # print "Ai_ui_P_vector(self.n0, q2theta_i(q, self.node_body_id)) =", Ai_ui_P_vector(self.n0, q2theta_i(q, self.node_body_id))
        #
        theta_edge = q2theta_i(q, self.node_body_id)
        
        self.r_iP = q2R_i(q, self.node_body_id) + Ai_ui_P_vector(self.r_iP[0:2], theta_edge)
        # print "n0 =", n0
        self.r_jP = q2R_i(q, self.edge_body_id) + Ai_ui_P_vector(self.r_jP[0:2], theta_edge)
        
        r_jPn0 = self.r_jP - self.r_jP
        # r_jRn0 = self.n0 - self.r_jR
        #
        # #   angle
        # fi_012 = self.angle(r_jPn0, -self.r_jRr_jP)
        # fi_021 = self.angle(r_jRn0, self.r_jRr_jP)
        r_jRr_jP_e = Ai_ui_P_vector(self.r_jRr_jP_e[0:2], theta_edge)
        
        #    
        self._distance_vector = (np.dot(r_jPn0, r_jRr_jP_e) * r_jRr_jP_e) - r_jPn0

        #   distance
        self._distance = np.linalg.norm(self._distance_vector, ord=2)

        self.normal_plus_distance = self.get_normal_2D() - self._distance_vector

        #   check if free node is inside (TF status)
        self._TF1 = np.linalg.norm(self.normal_plus_distance) < np.linalg.norm(self.normal)
        self._TF2 = self._distance <= self.max_penetration_depth
#         self._TF3 = fi_012 <= np.pi/2
#         self._TF4 = fi_021 <= np.pi/2
        
        #    distance sign
        if np.sign(np.cross(self.edge, r_jPn0)[2]) < 0:
            self._inside = False
            self._distance_sign = self._distance
        else:
            self._inside = True
            self._distance_sign = -self._distance

        #   assigne distance sign
        #   negative - indentation (contact)
        #   positive - separation (no contact)
#         if (np.linalg.norm(self.normal_plus_distance) < np.linalg.norm(self.normal)):# and (self._distance <= self.max_penetration_depth) and (fi_012 <= np.pi/2) and (fi_021 <= np.pi/2):
#             self._inside = True
#             self._distance_sign = -self._distance
#         else:
#             self._inside = False
#             self._distance_sign = abs(self._distance)

    def _vector_length(self, T1, T2):
        """
        Vector length is calculated as vector norm
        """
        #   vector formed between points T1, T2
        self._distance_vector = T2 - T1

        #   vector length
        _d = np.linalg.norm(self._distance_vector, ord=2)

        #   calculate normal for revolute clearance joint
        self.normal = self._distance_vector/_d
        # print "self.normal =", self.normal
        self.tanget = None

        #    check if node is inside or outside of an edge with normal n
        if np.sign(np.dot(self._distance_vector, self.normal)) <= -1:
            _inside = False
        else:
            _inside = True

        return _d, _inside

    def angle(self, vec1, vec2):
        return np.arccos(np.dot(vec1, vec2) / (np.linalg.norm(vec1, ord=2) * np.linalg.norm(vec2, ord=2)))

    def get_normal_2D(self):
        return self.normal[0:2]

    def get_tangent_2D(self):
        return self.tangent[0:2]

    def get_edge_2D(self):
        return np.array([self.edge_node_1[0:2], self.edge_node_2[0:2]])

    def get_node_2D(self):
        return np.array(self.r_iP[0:2])

    def distance_2D_in_direction_of_normal(self):
        """
        
        """
        self._distance_2D = 0

    def set_data_LCS(self, n0_LCS, r_jP_LCS, r_jR_LCS):
        """
        Function calculates contact data in LCS of each body in contact pair (node, edge)

        :return:
        """

    def contact_point_on_line(self):
        """
        Only implemented in subclass
        :return:
        """
        return None

if __name__ == "__main__":
    n0 = np.array([.08, .10, 0])
    r_jP = np.array([0.1, 0.1, 0])
    r_jR = np.array([.1, 0, 0])
    normal = np.array([1, 0, 0])
    # n0 = np.array([.08, .1, 0])
    # r_jP = np.array([0.1, 0.05, 0])
    # r_jR = np.array([.1, 0.1, 0])

    d = Distance(n0, r_jP, r_jR, normal)
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
    ax.arrow(r_jP[0], r_jP[1], r_jR[0]-r_jP[0], r_jR[1]-r_jP[1], width = 1E-4, head_width=4E-4, head_length=0.01, fc='k', ec='k')
    ax.plot([r_jP[0], r_jR[0]], [r_jP[1], r_jR[1]], color="black", marker="s", linestyle = "-")
    ax.text(r_jP[0]*1.01, r_jP[1]*1.01, "r_jP =["+str(r_jP[0])+","+str(r_jP[1])+"]", color='black')
    ax.text(r_jR[0]*1.01, r_jR[1]*1.01, "r_jR =["+str(r_jR[0])+","+str(r_jR[1])+"]", color='black')
    #   normal
    ax.arrow(r_jP[0], r_jP[1], normal[0]*1E-2, normal[1]*1E-2, width = 1E-4, head_width=4E-4, head_length=1E-4, fc='grey', ec='grey')

    #   plot n0r_jP
    ax.plot([n0[0], r_jP[0]], [n0[1], r_jP[1]], color="grey", linestyle = "--")

    plt.show()
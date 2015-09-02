'''
Created on 19. mar. 2015

@author: lskrinjar
'''
from pprint import pprint
import numpy as np

from ...q2q_body import q2q_body
from ...transform_cs import gcs2cm_lcs

class Distance(object):
    """
    classdocs
    """
    def __init__(self, n0, n1, n2=[], normal=[], node_body=None, edge_body=None, parent = None):
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

        if n2 == []:
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
        n1n0 = n1 - n0
        self.__tangent = n2n1 = n2 - n1
        self.tangent = n2n1_e = n2n1 / np.linalg.norm(n2n1)

        
        _distance_vector = (np.dot(n1n0, n2n1_e) * n2n1_e) - n1n0

        self._distance = np.linalg.norm(_distance_vector, ord=2)

        normal_plus_distance = self.normal + _distance_vector
        
        
        if np.linalg.norm(normal_plus_distance) < np.linalg.norm(self.normal):
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
    n0 = np.array([.5, +.1, 0])
    n1 = np.array([0, 0, 0])
    n2 = np.array([1, 0, 0])
    normal = np.array([0, 1, 0])
    d = Distance(n0, n1, n2, normal)
    pprint(vars(d))

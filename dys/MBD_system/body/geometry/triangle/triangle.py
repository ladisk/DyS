'''
Created on 7. feb. 2015

@author: lskrinjar
'''
import time
import numpy as np
from pprint import pprint
from MBD_system.body.geometry.edge.edge import Edge

class Triangle(object):
    '''
    classdocs
    '''
    def __init__(self, _normal, _id=0):
        '''
        Constructor
        '''
        #   triagnle id
        self._id = _id
        #   list of objects of type node
        self.nodes = []
        #   normal
        self.normal = np.array(_normal)
        
        
    def add_node(self, node):
        self.nodes.append(node)
    
    
    def intersects_plane(self, offset_z):
        """
        Returns True or False value if triangle intersects plane in z direction
        """
        cond_ = np.array([node.node[2] > offset_z for node in self.nodes]).any()  # .any()
        return cond_
    
    
    def plane_intersects_triangle(self, z_dim_lcs):
        """
        Method checks if triangle intersects plane parallel to z plane with offset z_dim_lcs
        and if plane intersects triangle it calculates intersection points - list self._edge
        and return true, if not it returns false
        :return:
        """
        #   predefined empty list of edge points
        self._edge = []
        for node1, node2 in zip(self.nodes, np.roll(self.nodes, 2, 0)):
            if min(node1.node[2], node2.node[2]) <= z_dim_lcs <= max(node1.node[2], node2.node[2]):

                # print node1.node
                t = (z_dim_lcs - node2.node[2]) / (node1.node[2] - node2.node[2])
                #   coordinates of intersection point
                P = node2.node+(node1.node - node2.node) * t

                #   append to list
                self._edge.append(P)

        if len(self._edge) == 2:
            return True
        else:
            return False


    def edge_in_area(self, _min, _max):
        """
        Function checks if list of line nodes (created from triangle plane intersection)
        :return:
        """
        #   predefined empy list to store booleans
        TF_list = []
        for _E in self._edge:
            for _E_xy, _min_xy, _max_xy in zip(_E, _min, _max):
                TF_list.append(_E_xy >= _min_xy)
                TF_list.append(_E_xy <= _max_xy)


        if all(TF_list) == True:
            return True
        else:
            return False



if __name__ == "__main__":
    a = Triangle(_normal=np.array([0, 1, 0]), _id=0)
    
    a.add_node(np.array([1., 0., 0.]))
#     a.add_node(np.array([1, 0, 2]))
    a.add_node(np.array([1., 0., 0.]))
    a.add_node(np.array([1., 0., 0.]))
#     pprint(vars(a))
    
#     a.intersects_plane(offset_z=1)

'''
Created on 15. feb. 2015

@author: lskrinjar
'''
import numpy as np
from pprint import pprint
class Node(object):
    '''
    classdocs
    '''
    def __init__(self, id_node, node, normal, id_triangle=None, parent=None):
        '''
        Constructor
        '''
        self._id = id_node
        self.parent = parent
        self.node = np.array(node)*1E-3
        self.normal = np.array(normal)
        self._id_triagle = id_triangle
    
    
    def to_2D(self):
        """
        
        """
        self.node_2D = self.node[0:2]
        self.normal_2D = self.normal[0:2]


if __name__ == "__main__":
    a = Node(id_node=0, node=np.array([1, 2, 3]), normal=np.array([.5, .5, 0]), id_triangle=0)
    a.node_2D()
    pprint(vars(a))

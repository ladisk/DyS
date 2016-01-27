'''
Created on 23. feb. 2015

@author: lskrinjar
'''
import numpy as np

class Edge(object):
    '''
    classdocs
    '''
    def __init__(self, node_i, node_j, normal):
        '''
        Constructor
        '''
        self.node_i = node_i
        self.node_j = node_j
        self.nodes = [self.node_i, self.node_j]
        self.normal_2D = normal
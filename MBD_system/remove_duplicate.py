'''
Created on 9. mar. 2015

@author: lskrinjar

'''
import numpy as np


def remove_duplicate(nodes):
    """
    Remove duplicate node vectors from matrix of nodes
    Args:
        nodes - a 2d array (or matrix) of nodes
    Raises:
        None
    Returns:
        unique_nodes - a 2d array (or matrix) of nodes without duplicates
        
    """
    nodes = np.ascontiguousarray(nodes)
    __unique_nodes = np.unique(nodes.view([('', nodes.dtype)] * nodes.shape[1]))
    unique_nodes = __unique_nodes.view(nodes.dtype).reshape((__unique_nodes.shape[0], nodes.shape[1]))
    
    return unique_nodes



'''
Created on 18. maj 2015

@author: luka.skrinjar
'''
import numpy as np
from d_ds import d_ds

def d2_ds2(node):
    """
    Function numerically calculates a second derivative d2/ds2 of node that defines a vector (line) 
    to calculate derivative
    """
    _d_ds = d_ds(node)
    
    
    _d2_ds2 = d_ds(_d_ds)
    
    
    return _d2_ds2


if __name__ == "__main__":
    _node = np.array([0.5, 0.2])
    
    print "_node =", _node
    print "d_ds =", d2_ds2(_node)
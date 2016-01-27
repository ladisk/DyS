'''
Created on 25. apr. 2015

@author: lskrinjar

'''
import numpy as np

def d_ds(node):
    """
    Function numerically calculates a derivative d/ds of a node that
    defines a vector to derive
    """
    if node[0] == 0:
        k = node[1]
    else:
        k = node[1] / node[0]
    
    _d_ds = np.array([1, k])

    return _d_ds


if __name__ == "__main__":
    _node = np.array([0.5, 0.2])
    
    print "_node =", _node
    print "d_ds =", d_ds(_node)
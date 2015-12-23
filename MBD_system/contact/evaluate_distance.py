__author__ = 'lskrinjar'

import numpy as np

def evaluate_distance_2D(n0, n1, n, t):
    """
    Function calculates a distance between a node and edge and returns
    a distance and inside/outside status if it is inside or outside of an edge
        n0 - free node in GCS
        n1 - edge node (start) in GCS
        normal - unit vector
        tanget - unit vector
    :return:
    """
    #   construct vector from free point to start of edge
    n1n0 = n1 - n0

    #    distance vector 
    distance_vector = (np.dot(n1n0, t) * t) - n1n0

    #    distance value
    d = np.linalg.norm(distance_vector, ord=2)

    #    check if node is inside or outside of an edge with normal n
    if np.sign(np.dot(n1n0, n)) <= -1:
        inside = False
    else:
        inside = True

    return d, inside
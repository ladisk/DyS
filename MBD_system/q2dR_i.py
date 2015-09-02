'''
Created on 30. maj 2014

@author: lskrinjar (email: skrinjar.luka@gmail.com)
'''
import numpy as np


def q2dR_i(q, body_id, N_b = None):
    """
    Function returns the velocity of a body from vector q and body_id.
    Args:
        q - vector of positions and velocities
        body_id - body_id generated when creating body object
        N_b - number of bodies
    Returns:
        dR_i_ - translational velocity of a body with id body_id
    Raises:
        None
    """
    if body_id == "ground":
        dR_i_ = 0
    else:
        if N_b is None:
            N_b = len(q) / 6
        dR_i_ = q[(3 * N_b) + 3 * body_id:(3 * N_b) + 3 * body_id + 2]
        
    return dR_i_

if __name__ == "__main__":
    N_b = 2
    print "N_b =", N_b
    q = np.arange(0, 6 * N_b)
    print "q =", q
#     print "q[8] =", q[8]
#     print "q[11] =", q[11]
    body_id = 0
    print "body_id =", body_id
    print q2dR_i(q, body_id, N_b)

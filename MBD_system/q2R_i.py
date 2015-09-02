'''
Created on 30. maj 2014

@author: lskrinjar (email: skrinjar.luka@gmail.com)
'''
import numpy as np


def q2R_i(q, body_id):
    """
    Function returns the angular velocity of a body from vector q and body_id.
    Args:
        q - vector of positions and velocities
        body_id - body_id generated when creating body object
        N_b - number of bodies
    Returns:
        dtheta_i_ - angular velocity of a body with id body_id
    Raises:
        None
    """
    if body_id == "ground":
        R_i = np.zeros(2)
    else:
        R_i = q[3 * body_id:3 * body_id + 2]
        
    return R_i

if __name__ == "__main__":
    N_b = 2
    print "N_b =", N_b
    q = np.arange(0, 6 * N_b)
    print "q =", q
    print "q[8] =", q[8]
    print "q[11] =", q[11]
    body_id = 1
    print q2R_i(q, body_id)

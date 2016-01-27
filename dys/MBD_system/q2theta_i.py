'''
Created on 21. apr. 2014

@author: lskrinjar (email: skrinjar.luka@gmail.com)
'''
import numpy as np


def q2theta_i(q, body_id):
    """
    Function returns the angle of rotation of a body from vector q and body_id.
    Args:
        q - vector of positions and velocities
        body_id - body_id generated when creating body object

    Returns:
        theta_i_ - angle of rotation of a body with id body_id
    Raises:
        None
    """
    if body_id == "ground":
        dtheta_i_ = 0
    else:
        dtheta_i_ = q[3*body_id+2]
        
    return dtheta_i_

if __name__ == "__main__":
    N_b = 2
    print "N_b =", N_b
    q = np.arange(0, 6*N_b)
    print "q =", q
    print "q[8] =", q[8]
    print "q[11] =", q[11]
    body_id = 0
    print q2theta_i(q, body_id)
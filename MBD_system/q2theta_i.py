"""
Created on 21. apr. 2014

@author: lskrinjar (email: skrinjar.luka@gmail.com)
"""
import numpy as np
from global_variables import GlobalVariables
from q2e_x_jk import q2e_x_jk


def q2theta_i(q, body_id, node_id=None):
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
        theta = 0.

    else:
        if hasattr(GlobalVariables, "q_i_dim"):
            #   rigid body
            if GlobalVariables.q_i_dim[body_id] == 3:
                theta = q[int(np.sum(GlobalVariables.q_i_dim[0:body_id])) + 2]

            #   point mass
            elif GlobalVariables.q_i_dim[body_id] == 2:
                theta = 0.

            #   flexible body
            else:
                e_x = q2e_x_jk(q, body_id=body_id, node_id=node_id)
                theta = np.arctan2(e_x[1], e_x[0])
        else:
            theta = 0.
        
    return theta

if __name__ == "__main__":
    N_b = 2
    print "N_b =", N_b
    q = np.arange(0, 6*N_b)
    print "q =", q
    print "q[8] =", q[8]
    print "q[11] =", q[11]
    body_id = 0
    print q2theta_i(q, body_id)
"""
Created on 21. apr. 2014

@author: lskrinjar (email: skrinjar.luka@gmail.com)
"""
import numpy as np
from global_variables import GlobalVariables
from q2de_x_jk import q2de_x_jk
from q2e_x_jk import q2e_x_jk
from ancf.ancf_omega import ancf_omega


def q2dtheta_i(q, body_id, node_id=None):
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
    #   ground
    if not isinstance(body_id, int):
        dtheta = 0.

    #   body
    else:
        if hasattr(GlobalVariables, "q_i_dim"):
            #   rigid body
            if GlobalVariables.q_i_dim[body_id] == 3:
                dtheta = q[int((len(q)/2) + np.sum(GlobalVariables.q_i_dim[0:body_id])) + 2]

            #   point mass
            elif GlobalVariables.q_i_dim[body_id] == 2:
                dtheta = 0.

            #   flexible body
            else:
                e_x = q2e_x_jk(q, body_id=body_id, node_id=node_id)
                de_x = q2de_x_jk(q, body_id=body_id, node_id=node_id)
                dtheta = ancf_omega(e_x, de_x)
        else:
            dtheta = 0.
    
    return dtheta


if __name__ == "__main__":
    N_b = 2
    print "N_b =", N_b
    q = np.arange(0, 6 * N_b)
    print "q =", q
    print "q[8] =", q[8]
    print "q[11] =", q[11]
    body_id = 1
    print "theta_i =", q2dtheta_i(q, body_id)

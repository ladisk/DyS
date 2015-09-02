'''
Created on 20. apr. 2015

@author: lskrinjar

'''
import numpy as np

from q2R_i import q2R_i
from q2theta_i import q2theta_i


def q2q_body(q, body_id):
    """
    Function returns the velocity vector of a body from vector q and body_id.
    Args:
        q - vector of positions and velocities
        body_id - body_id generated when creating body object
        N_b - number of bodies
    Returns:
        dR_i_ - velocity vector of translational and angular velocity of a body with id body_id
    Raises:
        None
    """
    q_body = np.zeros(3)

    q_body[0:2] = q2R_i(q, body_id)
    q_body[2] = q2theta_i(q, body_id)
    
    return q_body
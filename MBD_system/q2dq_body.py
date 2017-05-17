"""
Created on 20. apr. 2015

@author: lskrinjar

"""
import numpy as np


from global_variables import GlobalVariables


def q2dq_body(q, body_id, q_i_dim=None):
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
    if q_i_dim is None:
        q_i_dim = GlobalVariables.q_i_dim

    dq_body = q[np.sum(q_i_dim) + np.sum(q_i_dim[0:body_id]):np.sum(q_i_dim) + np.sum(q_i_dim[0:body_id + 1])]
    return dq_body


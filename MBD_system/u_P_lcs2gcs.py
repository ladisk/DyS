"""

created by: lskrinjar
date of creation: 25/01/2016
time of creation: 19:33
"""
from Ai_ui_P import Ai_ui_P_vector
from q2R_i import q2R_i
from q2theta_i import q2theta_i


def u_P_lcs2gcs(uP_lcs, q, body_id, node_id=None):
    """
    Function evaluates coordinates of a point uP in GCS based on coordinates in LCS
    :param uP_lcs:  a point vector in LCS 2D (x, y) (type: numpy array)
    :param q:       a vector of system coordinates and velocities (type:numpy array)
    :param id:      a id of a body (int)
    :return rP_gcs: a vector of point in GCS
    """
    rP_gcs = q2R_i(q, body_id) + Ai_ui_P_vector(uP_lcs, q2theta_i(q, body_id, node_id=node_id))
    return rP_gcs


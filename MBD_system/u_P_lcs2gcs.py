"""

created by: lskrinjar
date of creation: 25/01/2016
time of creation: 19:33
"""
from q2R_i import q2R_i
from q2theta_i import q2theta_i
from Ai_ui_P import Ai_ui_P_vector

def u_P_lcs2gcs(uP_lcs, q, id):
    """
    Function evaluates coordinates of a point uP in GCS based on coordinates in LCS
    :param uP_lcs:  a point vector in LCS 2D (x, y) (type: numpy array)
    :param q:       a vector of system coordinates and velocities (type:numpy array)
    :param id:      a id of a body (int)
    :return:
    """
    rP_gcs = q2R_i(q, id) + Ai_ui_P_vector(uP_lcs, q2theta_i(q, id))
    return rP_gcs

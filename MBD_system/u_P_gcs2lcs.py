"""

created by: lskrinjar
date of creation: 27/01/2016
time of creation: 17:49
"""
import numpy as np


from q2R_i import q2R_i
from q2theta_i import q2theta_i
from Ai_ui_P import Ai_ui_P_vector
from A_inv import A_inv_matrix


def u_P_gcs2lcs(rP_gcs, q, id):
    """
    Function evaluates coordinates of a point uP in LCS based on coordinates in GCS
    :param uP_lcs:  a point vector in LCS 2D (x, y) (type: numpy array)
    :param q:       a vector of system coordinates and velocities (type:numpy array)
    :param id:      a id of a body (int)
    :return:
    """
    uP_gcs = np.dot(A_inv_matrix(q2theta_i(q, id)), rP_gcs - q2R_i(q, id))
    return uP_gcs

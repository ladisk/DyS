'''
Created on 17. apr. 2014

@author: lskrinjar (email: skrinjar.luka@gmail.com)
'''
from cad2cm_lcs import cad2cm_lcs
import numpy as np


def u_P_cad2cm_lcs(_body_id, body, _u_P=np.zeros(2)):
    """
    Function assigns the body properties to body (from list of bodies) and calculates u_P vector in CM LCS of a body for each body in joint
    Args:
        bodies - list of bodies
    Returns:
        body_id_
        body_j_
        u_iP_f_cm_
        u_jP_f_cm_
    """
    #    check if id of body i is ground
    if _body_id == "ground":
        _u_P_f_cm = _u_P

    else:
        _u_P_f_cm = cad2cm_lcs(_u_P, body.CM_CAD_LCS[0:2])

    return _u_P_f_cm
'''
Created on 17. apr. 2014

@author: lskrinjar (email: skrinjar.luka@gmail.com)
'''
from cad2cm_lcs import cad2cm_lcs
import numpy as np


def ui_P_cad2cm_lcs(_body_id_i, _body_id_j, _u_iP, _u_jP, _bodies):
    """
    Function assigns the body properties to body (from list of bodies) and calculates u_P vector in CM LCS of a body for each body in joint
    Args:
        bodies - list of bodies
    Returns:
        body_i_
        body_j_
        u_iP_f_cm_
        u_jP_f_cm_
    """
    #    check if id of body i is ground
    if _body_id_i == "ground":
        _u_iP_f_cm = _u_iP
    else:
        _body_i = _bodies[int(_body_id_i)]
        _u_iP_f_cm = cad2cm_lcs(vector_in_CAD_LCS = _u_iP, CM_CAD_LCS = _body_i.CM_CAD_LCS[0:2])
    
    if _body_id_j == "ground":
        _u_jP_f_cm = _u_jP
    else:
        _body_j = _bodies[int(_body_id_j)]
        _u_jP_f_cm = cad2cm_lcs(vector_in_CAD_LCS = _u_jP, CM_CAD_LCS = _body_j.CM_CAD_LCS[0:2])
        
    return _u_iP_f_cm, _u_jP_f_cm
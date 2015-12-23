'''
Created on 20. apr. 2015

@author: lskrinjar

'''
import numpy as np

from A_inv import A_inv_matrix
from A import A_matrix
from cad2cm_lcs import cad2cm_lcs

def uP_gcs2lcs(u_P, theta):
    _u_P = np.dot(A_inv_matrix(theta), u_P)
    return _u_P

def uP_lcs2gcs(u_P, theta):
    _u_P = np.dot(A_matrix(theta), u_P)
    return _u_P

def gcs2cm_lcs(r_P, R, theta):
    """
    Function calculates a coordinate of a point P in GCS - global coordinate system
    from given coordinate in LCS - local coordinate system
    """
    if len(r_P) == 3:
        r_P = r_P[0:2]
    
    u_P = r_P - R
    
    _u_P = np.dot(A_inv_matrix(theta), u_P)
    return _u_P

def cm_lcs2gcs(_u_P, R, theta):
    """
    Function calculates a coordinate of a point P in LCS - local coordinate system
    from given coordinate in GCS - global coordinate system
    Args:
        r_P - position vector of point P in LCS
        R_i - 
        theta_i - 
    Returns:
        
    """
    if len(_u_P) == 3:
        _u_P = _u_P[0:2]
    u_P = np.dot(A_matrix(theta), _u_P)
    
    r_P = R + u_P
    return r_P

def gcs2lcs_z_axis(Rz_body, z_gcs):
    """
    Function tranforms z coordinates in global coordinate system to z coordinate in local coordinate system
    of a body
    :param Rz_body: a value of Z coordinate of body center of mass - local coordinate system
    :param z_gcs: a value of Z coordinate in GCS
    :return z_lcs: a coordinate of a point in LCS of a body
    """
    z_lcs = z_gcs - Rz_body
    return z_lcs

def u_P_cad2cm_lcs(_body_id, body, _u_P):
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
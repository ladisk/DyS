"""
Created on 20. apr. 2015

@author: lskrinjar

"""
import numpy as np

from A import A_matrix
from A_inv import A_inv_matrix
from cad2cm_lcs import cad2cm_lcs


def uP_gcs2lcs(r_P, theta):
    u_P = np.dot(A_inv_matrix(theta), r_P)
    return u_P


def uP_lcs2gcs(u_P, theta):
    u_P = np.dot(A_matrix(theta), u_P)
    return u_P


def gcs2cm_lcs(r_P, R=np.zeros(2), theta=0.):
    """
    Function calculates a coordinate of a point P in GCS - global coordinate system
    from given coordinate in LCS - local coordinate system
    """
    if len(r_P) == 3:
        r_P = r_P[0:2]
    
    u_P = r_P - R
    
    u_P = np.dot(A_inv_matrix(theta), u_P)
    return u_P


def cm_lcs2gcs(u_P, R, theta):
    """
    Function calculates a coordinate of a point P in LCS - local coordinate system
    from given coordinate in GCS - global coordinate system
    Args:
        r_P - position vector of point P in LCS
        R_i - 
        theta_i - 
    Returns:
        
    """
    if len(u_P) == 3:
        u_P = u_P[0:2]
    u_P = np.dot(A_matrix(theta), u_P)
    
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


def u_P_cad2cm_lcs(body_id, body, _u_P):
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
    if body_id == "ground":
        _u_P_f_cm = _u_P

    else:
        _u_P_f_cm = cad2cm_lcs(_u_P, body.u_CAD[0:2], 0)

    return _u_P_f_cm


if __name__ == "__main__":
    # r_P = np.array([-0.0042006637945, -0.002977900952])
    # theta = np.deg2rad(14.8)
    # R_i = np.array([-0.0020699993, -0.00227000055])
    # print "uPi =", gcs2cm_lcs(r_P, R=R_i, theta=theta)

    #   paper 2 flexible to rigid body data
    #   point of rigid constraint
    # r_P = np.array([-0.003482, 0.0034186])
    #   point of spring on flexible body
    r_P = np.array([-0.003482089781, 0.003481590506])
    R = np.array([-0.0042066, -0.0029779])
    theta = np.deg2rad(83.6)
    print "uPi =", gcs2cm_lcs(r_P, R=R, theta=theta)

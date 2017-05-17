"""
Created on 2. jun. 2015

@author: luka.skrinjar
"""
import numpy as np
from .A_theta import A_theta_matrix

def dr_contact_point_uP(dR, theta, dtheta, uP):
    """
    Function calculates velocity of contact point on a body
    Args:
        dR - translation velocity vector of a body (coordinate system)
        theta - body rotational position
        dtheta - body angular velocity
        uP - vector of contact point in body coordinate system (LCS)
    Returns:
        dr - velocity of (contact) point in GCS
    """
    # dr = dR + dtheta * (np.dot(A_theta_matrix(theta), uP))
    dr = dR + dtheta * np.dot(A_theta_matrix(theta), uP)
    return dr


if __name__ == "__main__":
    dR_OK = np.array([])
    theta_OK = 0.
    dtheta_OK = 0.
    uP_OK = np.array([])
    print dr_contact_point_uP(dR_OK, theta_OK, dtheta_OK, uP_OK)

    dR_NOK = np.array([])
    theta_NOK = 0.
    dtheta_NOK = 0.
    uP_NOK = np.array([])
    print dr_contact_point_uP(dR_NOK, theta_NOK, dtheta_NOK, uP_NOK)

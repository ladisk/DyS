'''
Created on 3. jun. 2014

@author: lskrinjar
'''

import itertools

import numpy as np

from MBD_system.A_theta import A_theta_matrix


class Force_Q_e_matrix(object):
    """
    classdocs
    """
    __id = itertools.cycle([0, 1])

    def __init__(self, body_id, u_iP, theta):
        """
        Create matrix that is used to create a vector of external force applied to a body
        Args:
            u_iP_f -
            theta -
        Returns:
            matrix_
        Raises:
            none
        """
        if body_id == -1 or body_id == "ground":
            self.matrix = np.zeros([3, 2])
        else:
            self.matrix = np.vstack((np.eye(2), np.dot(A_theta_matrix(theta), u_iP)))
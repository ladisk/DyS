"""

created by: lskrinjar
date of creation: 10/12/2016
time of creation: 22:21
"""
import numpy as np


from MBD_system.joint.joint import Joint
from MBD_system.joint.joint_C_q_matrix import Joint_C_q_matrix
from global_variables import GlobalVariables
from MBD_system.joint.support import Support
from global_variables import GlobalVariables


class SlopeDiscontinuity(Support):
    """
    classdocs
    """
    def __init__(self, body_id, node_id_i=0, node_id_j=0, properties_dict={}, parent=None):
        """

        :param body_id_i:
        :param body_id_j:
        :param u_iP_CAD:
        :param u_jP_CAD:
        :param u_iQ:
        :param parent:
        :return:
        """
        self.joint_type = "slope discontinuity"
        super(SlopeDiscontinuity, self).__init__(self.joint_type, body_id, body_id, node_id_i=node_id_i, node_id_j=node_id_j, properties_dict=properties_dict, parent=parent)

        self._parent = parent

        #   joint DOF
        self.joint_DOF = 1

        #   number of constrained nodal coordinates per node of finite element
        self.n_CNC = 1

        #   node id
        self.node_id = node_id_i

        #   body id
        self.body_id = body_id

        #   C_q matrix dimensions [rows, cols]
        # if self._parent._parent.bodies[body_id].q_i_size is None:
        #     self._parent._parent.bodies[body_id].evaluate_q_i_size()

        # self.C_q_dim = [self.n_CNC, self._parent._parent.bodies[body_id].q_i_size]
        self.C_q_dim = [self.n_CNC, 2]

    def evaluate_C_q(self, q):
        """

        :param q:
        :return:
        """
        self.C_q_list = []
        for body_id, node_id in zip(self.body_id_list, self.node_id_list):
            C_q_i = np.ones(2)

            self.C_q_list.append(C_q_i)

        return self.C_q_list

    def evaluate_Q_d(self, q):
        """

        :param q:
        :return:
        """
        Q_d = np.zeros(self.n_CNC)
        return Q_d

    def evaluate_C0(self, q):
        """

        :param q:
        :return:
        """


    def evaluate_C(self, q, t=None):
        """

        :param q:
        :param t:
        :return:
        """


    def evaluate_Q_c(self, L):
        """
        Function evaluates
        :param L:    vector of lagrange multipliers of a joint
        """
        return np.zeros(self.n_CNC)
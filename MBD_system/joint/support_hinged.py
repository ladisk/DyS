"""

created by: lskrinjar
date of creation: 13/08/2016
time of creation: 11:37
"""

# coding=utf-8


import numpy as np


from MBD_system.joint.joint_C_q_matrix import Joint_C_q_matrix
from MBD_system.joint.joint import Joint
from global_variables import GlobalVariables
from MBD_system.joint.support import Support
from MBD_system.u_P_lcs2gcs import u_P_lcs2gcs


class SupportHinged(Support):
    """
    classdocs
    """
    def __init__(self, body_id_i, body_id_j, node_id_i=0, node_id_j=0, properties_dict={}, parent=None):
        """

        :param joint_type:
        :param body_id_i:
        :param body_id_j:
        :param u_iP_CAD:
        :param u_jP_CAD:
        :param u_iQ:
        :param parent:
        :return:
        """
        #   type
        self.joint_type = "hinged support"
        super(SupportHinged, self).__init__(self.joint_type, body_id_i, body_id_j, node_id_i=node_id_i, node_id_j=node_id_j, properties_dict=properties_dict, parent=parent)

        #   joint DOF
        self.joint_DOF = 2

        #   number of constrained nodal coordinates per node of finite element
        self.n_CNC = 2

        #   size of vector q for body i and j
        self.e_n_i = GlobalVariables.q_i_dim[self.body_id_i]
        self.e_n_list = [self.e_n_i, None]
        self.C_q_dim = [self.n_CNC, self.e_n_list]

        #   initial evaluation of constraint equation for fixed joint
        # self.C0 = self.evaluate_C0(self.q0)
        # self.theta0 = self.C0[2]

    def evaluate_C_q(self, q=None):
        """
        Function evaluates C_q matrix for every element connected by this support (joint) object
        :param q:
        :return:
        """
        C_q_list = []
        #    construct C_q matrix for fixed support (joint)
        for body_id, node_id in zip(self.body_id_list, self.node_id_list):
            #   ground
            if type(body_id) is not int:
                C_q = Joint_C_q_matrix(self.joint_type, parent=self)

            #   actual body object
            else:
                C_q = Joint_C_q_matrix(self.joint_type, body_id=body_id, parent=self)

                C_q.matrix = self._parent._parent.bodies[body_id].evaluate_C_q_hinged(node_id)

            C_q_list.append(C_q)

        [C_qi, C_qj] = C_q_list

        return C_qi.matrix, C_qj.matrix

    def evaluate_Q_d(self, q):
        """

        :param q:
        :return:
        """
        Q_d = np.zeros(2)
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
        for i, (body_id, node_id, uP) in enumerate(zip(self.body_id_list, self.node_id_list, self.u_P_LCS_list)):
            #   if body
            if isinstance(body_id, int):
                rP = self._parent._parent.bodies[body_id].evaluate_r(q, node_id=node_id)
            #   if ground
            else:
                rP = uP

            self.r_P_GCS_list[i] = rP

        [rP_i, rP_j] = self.r_P_GCS_list

        C = rP_i - rP_j

        return C
"""

created by: lskrinjar
date of creation: 24/12/2016
time of creation: 11:53
"""
import numpy as np


from MBD_system.joint.joint import Joint
from MBD_system.joint.support import Support
from MBD_system.Ai_theta_ui_P import Ai_theta_ui_P_vector
from MBD_system.q2theta_i import q2theta_i
from MBD_system.ancf.mesh.finite_element.beam_2d_euler_bernoulli import Beam2DEulerBernoulli
from MBD_system.ancf.mesh.finite_element.beam_2d_shear_deformable import Beam2DShearDeformable
from MBD_system.Ai_ui_P import Ai_ui_P_vector
from MBD_system.q2dtheta_i import q2dtheta_i
from MBD_system.q2de_jk import q2de_jk
from MBD_system.u_P_lcs2gcs import u_P_lcs2gcs
from MBD_system.q2q_body import q2q_body
from global_variables import GlobalVariables


class JointRevoluteFlexibleFlexible(Support):
    """
    Class of revolute joint between a two flexible bodies
    Ref: Describing Rigid-Flexible Multibody Systems Using Absolute Coordinates (doi:10.1023/B:NODY.0000014553.98731.8d)
    """
    def __init__(self, body_id_i, body_id_j, node_id_i=None, node_id_j=None, properties_dict={}, parent=None):
        """
        :param body_id_i:   id of rigid body
        :param body_id_j:   id of flexible body
        :param node_id_i:
        :param node_id_j:
        :return:
        """
        self.joint_type = "revolute joint flexible-flexible"
        super(JointRevoluteFlexibleFlexible, self).__init__(self.joint_type, body_id_i, body_id_j, node_id_i=node_id_i, node_id_j=node_id_j, parent=parent)

        #   number of constrained nodal coordinates per node of finite element
        self.n_CNC = 2

        #   size of vector q for body i and j
        self.e_n_i = GlobalVariables.q_i_dim[self.body_id_i]
        self.e_n_j = GlobalVariables.q_i_dim[self.body_id_j]
        self.e_n_list = [self.e_n_i, self.e_n_j]

        #   C_q matrix dimensions [rows, cols]
        self.C_q_dim = [self.n_CNC, self.e_n_list]

    def evaluate_Q_d(self, q):
        """
        Function is defined in subclass
        :param q:
        :return:
        """
        # u_P_i = self.u_P_LCS_list[self.body_id_list.index(self.body_id_i)]
        # theta_i = q2theta_i(q, self.body_id_i)
        # Q_d_i = Ai_ui_P_vector(u_P_i, theta_i) * (q2dtheta_i(q, self.body_id_i) ** 2)
        #
        # Q_d_j = np.zeros(self.n_CNC)

        Q_d = np.zeros(self.n_CNC)
        return Q_d

    def evaluate_C0(self, q):
        """
        Function is defined in subclass
        :param q:
        :return:
        """
        return None

    def evaluate_C(self, q, t=None):
        """
        Function evaluates C for the joint
        :param q:
        :param t:
        :return:
        """
        for i, (body_id, node_id) in enumerate(zip(self.body_id_list, self.node_id_list)):
            self.r_P_GCS_list[i] = self._parent._parent.bodies[body_id].evaluate_r(q, node_id=node_id)

        [rP_i, rP_j] = self.r_P_GCS_list

        self.C = rP_i - rP_j

        return self.C

    def evaluate_C_q(self, q):
        """

        :param q:
        :return:
        """
        if not self.C_q_list:
            #   evaluate constraint matrix only once, as it is constant
            for body_id, node_id, sgn, n_e_i in zip(self.body_id_list, self.node_id_list, [+1, -1], self.C_q_dim[1]):
                C_q = self._evaluate_C_q(n_e_i, body_id, node_id, sgn)
                self.C_q_list.append(C_q)

        #   create list
        [C_q_i, C_q_j] = self.C_q_list

        return C_q_i, C_q_j

    def _evaluate_C_q(self, n_e, body_id, node_id, sgn):
        """
        Evaluate constraint jacobian matrix for flexible body
        :param n_e:     number of nodal coordinates of a body (length of vector e)
        :return:
        """
        C_q = np.zeros([2, n_e])

        node_e_n = self._parent._parent.bodies[body_id].mesh.node_e_n

        C_q[0:2, node_id * node_e_n:node_id * node_e_n + 2] = sgn * np.eye(2)

        return C_q
"""

created by: lskrinjar
date of creation: 14/12/2016
time of creation: 14:27
"""
from matplotlib import pyplot as plt
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
from MBD_system import transform_cs
from MBD_system.q2q_body import q2q_body
from MBD_system.q2R_i import q2R_i
from global_variables import GlobalVariables
from simulation_control_widget.vtk_widget.marker.marker import Marker
from MBD_system.extract_from_dictionary_by_string_in_key import extract_from_dictionary_by_string_in_key
from MBD_system.q2e_x_jk import q2e_x_jk
from MBD_system.q2de_x_jk import q2de_x_jk
from MBD_system.r_ij_P import r_ij_P


class JointRigidPointMassFlexible(Joint):
    """
    Class of rigid joint between a rigid body and a flexible body
    Ref: Describing Rigid-Flexible Multibody Systems Using Absolute Coordinates (doi:10.1023/B:NODY.0000014553.98731.8d)
    """

    def __init__(self, body_id_i, body_id_j, u_iP=np.zeros(2, dtype=float), node_id_j=None, properties_dict={}, parent=None):
        """
        :param support:
        :param body_id_i:   id of point mass
        :param body_id_j:   id of flexible body
        :param u_iP_CAD:
        :param u_jP_CAD:
        :param u_iQ:
        :param parent:
        :return:
        """
        self.joint_type = "rigid joint point mass-flexible"
        super(JointRigidPointMassFlexible, self).__init__(self.joint_type, body_id_i, body_id_j, u_iP_CAD=u_iP, properties_dict=properties_dict, parent=parent)

        #   number of constrained nodal coordinates per node of finite element
        self.n_CNC = 2

        #   body type list
        self.body_type_list = ["point mass", "flexible"]

        #   size of vector q for body i and j
        self.e_n_i = GlobalVariables.q_i_dim[self.body_id_i]
        self.e_n_j = GlobalVariables.q_i_dim[self.body_id_j]
        self.e_n_list = [self.e_n_i, self.e_n_j]

        #   C_q matrix dimensions [rows, cols]
        self.C_q_dim = [self.n_CNC, self.e_n_list]
        self.C_q_i = None
        self.C_q_j = None

        #   node id (node k of body j)
        self.node_id_j = node_id_j

        #   gradient vector on flexible body j
        self.e_j = self._parent._parent.bodies[self.body_id_j].e_node(self.q0, self.node_id_j)[0:2]
        self.e_j_grad = self._parent._parent.bodies[self.body_id_j].e_node(self.q0, self.node_id_j)[2:4]

        #   vector of joint position in point mass coordinate system

        if (u_iP == np.zeros(2)).all():
            u_iP = self.e_j - q2R_i(self.q0, self.body_id_i)

        self.u_iP_LCS = self.u_P_LCS_list[0] = self.u_iP = u_iP

        #   list of markers
        self.markers, self.r_P_GCS_list = self._create_markers()

        #   update lists
        self._update_lists()

    def _reset(self, q):
        """

        :return:
        """
        #   node on rigid body in GCS
        rP_i = u_P_lcs2gcs(self.u_P_LCS_list[self.body_id_list.index(self.body_id_i)], q, self.body_id_i)

        #   node on flexible body
        rP_j = self._parent._parent.bodies[self.body_id_j].evaluate_r(q, node_id=self.node_id_j)

        #   list of coordinates vectors
        self.r_P_GCS_list = [rP_i, rP_j]

    def evaluate_rijP(self, q):
        """

        :param q:
        :return:
        """
        rijP = (q2R_i(q, self.body_id_i) + self.u_iP_LCS) - self._parent._parent.bodies[self.body_id_j].e_node(q, self.node_id_j)[0:2]
        return np.linalg.norm(rijP, ord=2)

    def evaluate_Q_d(self, q):
        """
        Function is defined in subclass
        :param q:
        :return:
        """
        Q_d = np.zeros(self.n_CNC)

        # #   get uPi in LCS
        # _u_P_i = self.u_P_LCS_list[self.body_id_list.index(self.body_id_i)]
        #
        # #   rigid body i data
        # theta_i = q2theta_i(q, self.body_id_i)
        # omega_i = q2dtheta_i(q, self.body_id_i)
        #
        # #   evaluate uPi in GCS
        # u_P_i = Ai_ui_P_vector(_u_P_i, theta_i)
        # #   vector Qd
        # Q_d[0:2] = u_P_i * (omega_i ** 2)
        #
        # # Q_d[2] = (omega_i ** 2) * np.dot(u_P_i, self.e_j_grad)
        # e_x = q2e_x_jk(q, self.body_id_j, self.node_id_j)
        # Ai_uiR = Ai_ui_P_vector(self.u_iR, theta_i)
        # de_x = q2de_x_jk(q, self.body_id_j, self.node_id_j)
        # u_i = Ai_theta_ui_P_vector(self.u_iP_LCS, q2theta_i(q, self.body_id_i), 0)
        # Q_d[2] = np.dot(u_i, e_x) * omega_i + np.dot(u_P_i, de_x)

        # Q_d[2] = (np.dot(Ai_uiR, e_x) * (omega_i ** 2)) - (2. * omega_i * np.dot(Ai_theta_ui_P_vector(self.u_iR, theta_i, 0), de_x))

        return Q_d

    def evaluate_C0(self, q0, t=0.):
        """
        Function is defined in subclass
        :param q:
        :return:
        """
        self.C0 = self.evaluate_C(q0)

        return self.C0

    def evaluate_C(self, q, t=0.):
        """
        Function is defined in subclass
        :param q:
        :param t:
        :return:
        """
        self.C = np.zeros(self.n_CNC)
        #   node on rigid body in GCS
        rP_i = u_P_lcs2gcs(self.u_P_LCS_list[self.body_id_list.index(self.body_id_i)], q, self.body_id_i)

        #   node on flexible body
        rP_j = self._parent._parent.bodies[self.body_id_j].evaluate_r(q, node_id=self.node_id_j)

        #   list of coordinates vectors
        self.r_P_GCS_list = [rP_i, rP_j]

        #   vector C
        self.C[0:2] = rP_i - rP_j

        #   vector perpendicular to gradient vector
        # theta_i = q2theta_i(q, self.body_id_i)
        # rR_i = Ai_ui_P_vector(self.u_iR, theta_i)

        #   gradient vector at node of a joint
        # self.e_j_grad = self._parent._parent.bodies[self.body_id_j].e_node(q, self.node_id_j)[2:4]
        # self.e_j_grad = q2e_x_jk(q, self.body_id_j, self.node_id_j)
        # e_j_grad_unit = e_j_grad / np.linalg.norm(e_j_grad)
        #   scalar product of two vectors

        # rR_i_scaled = rR_i * np.linalg.norm(self.e_j_grad)
        # print "rR_i =", rR_i
        # print "self.e_j_grad =", self.e_j_grad
        # print "np.dot() =", np.dot(rR_i, self.e_j_grad)
        # print "fi =", np.rad2deg(np.arccos(np.dot(rR_i, self.e_j_grad) / (np.linalg.norm(rR_i) * np.linalg.norm(self.e_j_grad))))

        # plt.plot([0, rR_i[0]], [0, rR_i[1]], color="blue")
        # plt.plot([0, self.e_j_grad[0]], [0, self.e_j_grad[1]], color="green")
        # plt.show()

        # C[2] = np.dot(rR_i, e_j_grad)
        #   vector C
        # self.C[2] = np.dot(rR_i, self.e_j_grad)
        #   override for constraint check - this equations is not explicitly used at each integration time step
        # self.C[2] = 0.
        # C[2] = rR_i_scaled - self.e_j_grad

        # print "test =", np.dot(self.e_j_grad, Ai_ui_P_vector(self.u_iR, theta_i + np.deg2rad(90.)))

        # print "rR_i =", rR_i, "e_j_grad =", e_j_grad, "C[2:4] =", C[2:4], "rR_i - e_j_grad =", rR_i - e_j_grad, "dfi =", np.arccos(np.dot(rR_i, e_j_grad) / (np.linalg.norm(rR_i) * np.linalg.norm(e_j_grad))), np.linalg.norm(rR_i), np.linalg.norm(e_j_grad)

        # print "dfi =", np.rad2deg(np.arccos(np.dot(rR_i, e_j_grad_unit) / (np.linalg.norm(rR_i) * np.linalg.norm(e_j_grad_unit))))
        return self.C

    def evaluate_C_q(self, q):
        """

        :param q:
        :return:
        """
        #   gradient vector at node of a joint
        self.e_j_grad = self._parent._parent.bodies[self.body_id_j].e_node(q, self.node_id_j)[2:4]

        #   constraint jacobian matrix of a rigid body
        C_q_i = self._evaluate_C_q_i(q)

        #   constraint jacobian matrix for flexible body
        if self.C_q_j is None:
            self.C_q_j = self._evaluate_C_q_j(q)

        #   create list
        self.C_q_list = [C_q_i, self.C_q_j]

        return C_q_i, self.C_q_j

    def _evaluate_C_q_i(self, q):
        """
        Jacobian matrix of constraint equation of point mass body i
        :param q:
        :return:
        """
        C_q_i = np.zeros([self.C_q_dim[0], self.C_q_dim[1][0]])

        # C_q_i[0:2, 2] = Ai_theta_ui_P_vector(self.u_iP_LCS, q2theta_i(q, self.body_id_i), 0)

        C_q_i[0:2, 0:2] = np.eye(2)

        # C_q_i[2:4, 2] = Ai_theta_ui_P_vector(self.u_iR * np.linalg.norm(self.e_j_grad), q2theta_i(q, self.body_id_i), 0)
        # C_q_i[2, 2] = np.dot(self.e_j_grad, C_q_i[0:2, 2])

        #   update of code
        # ex_j = q2e_x_jk(q, body_id=self.body_id_j, node_id=self.node_id_j)
        # Ai_theta_uiR = Ai_theta_ui_P_vector(self.u_iR, q2theta_i(q, self.body_id_i), 0)
        # C_q_i[2, 2] = np.dot(Ai_theta_uiR, ex_j)

        return C_q_i

    def _evaluate_C_q_j(self, q):
        """
        Jacobian matrix of constraint equation of flexible body j
        :param q:
        :return:
        """
        C_q_j = np.zeros([self.n_CNC, self.e_n_j])

        #   constrain of position
        node_e_n = self._parent._parent.bodies[self.body_id_j].mesh.node_e_n
        C_q_j[0:2, self.node_id_j * node_e_n:self.node_id_j * node_e_n + 2] = -np.eye(2)
        # C_q_j[0:2, self.node_id_j * node_e_n + 2:self.node_id_j * node_e_n + 4] = -np.eye(2)

        return C_q_j

    def _create_markers(self):
        """
        Function create markers
        :return:
        """
        markers = []
        r_P_GCS_list = []
        node_id = None
        for i, (body_id, body) in enumerate(zip(self.body_id_list, self.body_list)):
            #   rigid body
            if body_id == self.body_id_i:
                uP = self.u_P_LCS_list[0]
                theta0 = np.zeros(3)
                rP = body.R[0:2] + uP

            #   flexible body
            else:
                uP = None
                rP = body.evaluate_r(self.q0, node_id=self.node_id_j)
                node_id = self.node_id_j
                e_j_grad = body.mesh.e_j(self.node_id_j)[2:4]
                theta0 = np.zeros(3)
                theta0[2] = np.arctan2(e_j_grad[1], e_j_grad[0])

            r_P_GCS_list.append(rP)

            marker = Marker(rP, uP=uP, theta0=theta0, body_id=body_id, node_id=node_id, body=body, parent=self)
            markers.append(marker)

        return markers, r_P_GCS_list

    def update_vtk_data(self, q):
        """

        :param q:
        :return:
        """
        for marker, rP in zip(self.markers, self.r_P_GCS_list):
            marker.update_vtk_data(q, rP=rP)


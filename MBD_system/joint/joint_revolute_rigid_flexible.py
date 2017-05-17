"""

created by: lskrinjar
date of creation: 14/12/2016
time of creation: 14:27
"""
import numpy as np


from MBD_system.joint.joint import Joint
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
from simulation_control_widget.vtk_widget.marker.marker import Marker
from MBD_system.q2e_x_jk import q2e_x_jk
from MBD_system.q2de_x_jk import q2de_x_jk
from MBD_system.q2e_jk import q2e_jk
from MBD_system import transform_cs
from MBD_system.q2R_i import q2R_i


class JointRevoluteRigidFlexible(Joint):
    """
    Class of rigid joint between a rigid body and a flexible body
    Ref: Describing Rigid-Flexible Multibody Systems Using Absolute Coordinates (doi:10.1023/B:NODY.0000014553.98731.8d)
    """

    def __init__(self, body_id_i, body_id_j, u_iP=np.array([0, 0], dtype=float), node_id_j=None, properties_dict={}, parent=None):
        """
        :param support:
        :param body_id_i:   id of rigid body
        :param body_id_j:   id of flexible body
        :param u_iP_CAD:
        :param u_jP_CAD:
        :param u_iQ:
        :param parent:
        :return:
        """
        self.joint_type = "revolute joint rigid-flexible"
        super(JointRevoluteRigidFlexible, self).__init__(self.joint_type, body_id_i, body_id_j, u_iP_CAD=u_iP, properties_dict=properties_dict, parent=parent)

        #   number of constrained nodal coordinates per node of finite element
        self.n_CNC = 2

        #   body type list
        self.body_type_list = ["rigid", "flexible"]

        #   size of vector q for body i and j
        self.e_n_i = GlobalVariables.q_i_dim[self.body_id_i]
        self.e_n_j = GlobalVariables.q_i_dim[self.body_id_j]
        self.e_n_list = [3, self.e_n_j]

        #   C_q matrix dimensions [rows, cols]
        #   number of columns is different for rigid and flexible body
        self.C_q_dim = [self.n_CNC, self.e_n_list]
        self.C_q_i = None
        self.C_q_j = None

        #   node id (node k of body j)
        self.node_id_j = node_id_j

        self.e_j_grad = q2e_x_jk(self.q0, self.body_id_j, self.node_id_j)
        self.e_j = q2e_jk(self.q0, self.body_id_j, self.node_id_j)

        #   list of markers
        self.markers = self._create_markers()

        #   update lists
        self._update_lists()

        #   number of nodal coordinates of a flexible body (mesh)
        self.e_n = self._parent._parent.bodies[self.body_id_j].mesh.n_NC

        #   vector perpendicular to gradient vector of ANCF flexible body j
        if (u_iP == np.zeros(2)).all():
            #   vector perpendicular to gradient vector
            r_iP = self.e_j

            R_i = q2R_i(self.q0, self.body_id_i)
            theta_i = q2theta_i(self.q0, self.body_id_i)
            self.u_iP = transform_cs.gcs2cm_lcs(r_iP, R=R_i, theta=theta_i)

        else:
            self.u_iP = u_iP

    def _reset(self, q):
        """

        :param q:
        :return:
        """
        #   node on rigid body in GCS
        rP_i = u_P_lcs2gcs(self.u_P_LCS_list[self.body_id_list.index(self.body_id_i)], q, self.body_id_i)

        #   node on flexible body
        rP_j = self._parent._parent.bodies[self.body_id_j].evaluate_r(q, node_id=self.node_id_j)

        #   list of coordinates vectors
        self.r_P_GCS_list = [rP_i, rP_j]

    def _create_markers(self):
        """
        Function create markers
        :return:
        """
        markers = []
        for i, (body_id, rP, body, node_id) in enumerate(zip(self.body_id_list, self.r_P_GCS_list, self.body_list, [self.node_id_i, self.node_id_j])):
            #   rigid body
            if body_id == self.body_id_i:
                uP = self.u_P_LCS_list[0]
                theta0 = np.zeros(3)
                theta0[2] = np.arctan2(self.e_j_grad[1], self.e_j_grad[0]) - body.theta[2]
            #   flexible body
            else:
                uP = None
                rP = body.evaluate_r(self.q0, node_id=self.node_id_j)

            marker = Marker(rP, uP=uP, theta0=theta0, body_id=body_id, node_id=node_id, body=body, parent=self)
            markers.append(marker)

        return markers

    def evaluate_Q_d(self, q):
        """
        Function is defined in subclass
        :param q:
        :return:
        """
        Q_d = np.zeros(self.n_CNC)

        u_P_i = self.u_P_LCS_list[self.body_id_list.index(self.body_id_i)]
        theta_i = q2theta_i(q, self.body_id_i)
        Q_d_i = Ai_ui_P_vector(u_P_i, theta_i) * (q2dtheta_i(q, self.body_id_i) ** 2)

        Q_d_j = np.zeros(self.n_CNC)

        Q_d = Q_d_i - Q_d_j
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
        #   node on rigid body in GCS
        rP_i = u_P_lcs2gcs(self.u_P_LCS_list[self.body_id_list.index(self.body_id_i)], q, self.body_id_i)

        #   node on flexible body
        rP_j = self._parent._parent.bodies[self.body_id_j].evaluate_r(q, node_id=self.node_id_j)

        #   list of coordinates vectors
        self.r_P_GCS_list = [rP_i, rP_j]

        #   vector C
        self.C = rP_i - rP_j

        return self.C

    def evaluate_C_q(self, q):
        """

        :param q:
        :return:
        """
        #   constraint jacobian matrix for rigid body
        C_q_i = self._evaluate_C_q_i(q)

        #   constraint jacobian matrix for flexible body
        C_q_j = self._evaluate_C_q_j(q)

        #   create list
        self.C_q_list = [C_q_i, C_q_j]

        return C_q_i, C_q_j

    def _evaluate_C_q_i(self, q):
        """
        Evaluate constraint jacobian matrix for rigid body i
        :param q:
        :return:
        """
        C_q_i = np.zeros([self.n_CNC, 3])
        C_q_i[0:2, 0:2] = np.eye(2)
        u_P_i = self.u_P_LCS_list[self.body_id_list.index(self.body_id_i)]

        theta_i = q2theta_i(q, self.body_id_i)

        C_q_i[:, -1] = Ai_theta_ui_P_vector(u_P_i, theta_i, self.body_id_list.index(self.body_id_i))

        return C_q_i

    def _evaluate_C_q_j(self, q):
        """
        Evaluate constraint jacobian matrix for flexible body j
        :param q:
        :return:
        """
        C_q_j = np.zeros([self.n_CNC, self.e_n_j])

        # if self._parent._parent.bodies[self.body_id_j].mesh.element_type == Beam2DEulerBernoulli.get_element_type():
        #     C_q_j = np.zeros([2, self.e_n])
        # if self._parent._parent.bodies[self.body_id_j].mesh.element_type == Beam2DShearDeformable.get_element_type():
        #     C_q_j = np.zeros([4, self.e_n])
        #
        # # print "Sj ="
        # #   vector of absolute nodal coordinates of a flexible body (of a mesh)
        # # e_j = q2q_body(q, self.body_id_j)
        # # print self._parent._parent.bodies[self.body_id_j].mesh.evaluate_S(e_j, self.node_id_jk)
        #
        # C_q_j[0:2, 0:2] = -np.eye(2)

        #   constrain of position
        node_e_n = self._parent._parent.bodies[self.body_id_j].mesh.node_e_n
        C_q_j[0:2, self.node_id_j * node_e_n:self.node_id_j * node_e_n + 2] = -np.eye(2)

        return C_q_j

    def update_vtk_data(self, q):
        """

        :param q:
        :return:
        """
        for marker, rP in zip(self.markers, self.r_P_GCS_list):
            marker.update_vtk_data(q, rP=rP)


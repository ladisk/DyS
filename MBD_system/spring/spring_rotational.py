"""
Created on 21. feb. 2014

@author: lskrinjar
"""
import itertools
from pprint import pprint
import numpy as np

from MBD_system.Ai_ui_P import Ai_ui_P_vector
from MBD_system.MBD_system_items import SpringItem
from MBD_system.dr_ij_P_dq import dr_ij_P_dq
from MBD_system.force.force import Force
from MBD_system.force.force_matrix import Force_Q_e_matrix
from MBD_system.force.force_vector import Force_Q_e_vector
from MBD_system.q2R_i import q2R_i
from MBD_system.q2dR_i import q2dR_i
from MBD_system.q2dtheta_i import q2dtheta_i
from MBD_system.q2theta_i import q2theta_i
from MBD_system.r_ij_P import r_ij_P
from MBD_system.u_P_cad2cm_lcs import u_P_cad2cm_lcs
from MBD_system.u_P_lcs2gcs import u_P_lcs2gcs
from simulation_control_widget.vtk_widget.marker.marker import Marker
from MBD_system.spring.spring import Spring
from global_variables import GlobalVariables
from MBD_system.q2e_jk import q2e_jk


class SpringRotational(Spring):
    """
    classdocs
    """
    def __init__(self, _name, body_id_i, body_id_j, u_iP_CAD=np.zeros(2), u_jP_CAD=np.zeros(2), u_iP_LCS=np.zeros(2, dtype=float), u_jP_LCS=np.zeros(2, dtype=float), k_t=0., c_t=0., theta_0=0., joint_id=None, properties_dict={}, parent=None):
        """

        :return:
        """
        super(SpringRotational, self).__init__(_name, body_id_i, body_id_j, u_iP_CAD=u_iP_CAD, u_jP_CAD=u_jP_CAD, u_iP_LCS=u_iP_LCS, u_jP_LCS=u_jP_LCS, properties_dict=properties_dict, joint_id=joint_id, parent=parent)

        #    spring type
        self.spring_type = "rotational"

        #   body type list
        self.body_type_list = ["rigid", "rigid"]

        #   add additional properties
        self.add_attributes_from_dict(self.properties)

        #   undeformed angle
        if hasattr(self._parent, "q0"):
            self.q0 = self._parent.q0

        #   add attributes from parent object
        self._add_attributes_from_parent()

        #   spring properties
        self.k_t = k_t
        self.c_t = c_t

        self.theta_0 = self._evaluate_theta(self.q0)

        #   time dependent variables
        self.theta = self.theta_0
        self.dtheta = 0.
        self.M_s = 0.

        #   create list of force objects
        self.Q_e_list = self._create_Q_e()

    def _evaluate_theta(self, q):
        """

        :return:
        """
        theta_list = np.zeros(2)
        for i, (body_id, node_id) in enumerate(zip(self.body_id_list, self.node_id_list)):
            if body_id == "ground":
                theta_body = 0.

            else:
                theta_body = q2theta_i(q, body_id, node_id=node_id)

            #    append theta of a body to array of both body theta angles
            theta_list[i] = theta_body
        # print "theta_list =", theta_list
        theta = theta_list[0] - theta_list[1]

        return theta

    def _evaluate_dtheta(self, q):
        """

        :param q:
        :return:
        """
        dtheta_list = np.zeros(2)
        for i, (body_id, node_id) in enumerate(zip(self.body_id_list, self.node_id_list)):
            if body_id == "ground":
                dtheta_body = 0.

            else:
                dtheta_body = q2dtheta_i(q, body_id, node_id=node_id)

            #    append theta of a body to array of both body angular velocities
            dtheta_list[i] = dtheta_body

        dtheta = dtheta_list[0] - dtheta_list[1]

        return dtheta

    def _evaluate_M_s(self, theta, dtheta):
        """

        :return:
        """
        return self.k_t * (theta - self.theta_0) + self.c_t * dtheta

    def evaluate_rijP(self, q):
        """

        :return:
        """
        r_P_list = []
        for body_id, uP, node_id in zip(self.body_id_list, self.u_P_LCS_list, self.node_id_list):
            if GlobalVariables.q_i_dim[body_id] == 3:
                rP = u_P_lcs2gcs(uP, q, body_id)

            else:
                # rP = q2R_i(q, body_id)
                rP = q2e_jk(q, body_id, node_id)

            r_P_list.append(rP)

        rijP = r_P_list[0] - r_P_list[1]

        return np.linalg.norm(rijP)

    def evaluate_Q_e(self, q):
        """
        Function evaluates a generalized external force vector
        :param q:
        :return:
        """
        #   theta
        self.theta = self._evaluate_theta(q)
        # print "self.theta =", self.theta
        #   dtheta
        self.dtheta = self._evaluate_dtheta(q)

        #   torsional force-moment
        self.M_s = self._evaluate_M_s(self.theta, self.dtheta)

        # Q_e_list = []
        for i, (body_id, Q_e) in enumerate(zip(self.body_id_list, self.Q_e_list)):
            if body_id != "ground":
                #   size of vector Q_e
                n_Q_e = GlobalVariables.q_i_dim[body_id]

                #   rigid (planar) body
                if n_Q_e == 3:
                    Q_e.Q_e[-1] = self.signs[self.body_id_list.index(int(body_id))] * self.M_s
                    Q_e.evaluate_Q_e(M=self.M_s)

                #   flexible (planar) body
                else:
                    # name = "spring_force_on_flexible_body_of_" + self._parent._name
                    # Q_e = Force(body_id, force_name=name, node_id=self._parent.node_id_list[i], element_id=self._parent.element_id_list[i], element_ksi=self._parent.element_ksi_list[i], parent=self)
                    Q_e.Mz = self.signs[self.body_id_list.index(int(body_id))] * self.M_s
                    Q_e.Q_e = Q_e._evaluate_Q_e_flexible_M(0., q, Q_e.Mz)

                    # Q_e.evaluate_Q_e(0., q)
            else:
                Q_e = Force_Q_e_vector()

            #    append spring force vector object to list
            # Q_e_list.append(Q_e)

        [Q_i, Q_j] = self.Q_e_list

        #   update force values in object (for data tracking and visualization)
        self._update_F(Q_e_list=self.Q_e_list)

        return Q_i.Q_e, Q_j.Q_e

    def _track_data(self, step, h, t):
        """

        :param step:
        :param h:
        :param t:
        :return:
        """

    def mechanical_energy(self, q):

        """

        :param q:
        :return:
        """
        return 0

    def print_F(self):
        """
        Function defined in subclass
        :return:
        """
        print "M ="
        print self._evaluate_M_s(self.theta, self.dtheta)
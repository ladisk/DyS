__author__ = 'lskrinjar'

import time
import numpy as np

from MBD_system.contact.contact import Contact
from MBD_system.q2R_i import q2R_i
from MBD_system.contact.distance.distance import Distance
from MBD_system.A import A_matrix
from MBD_system.dr_contact_point_uP import dr_contact_point_uP
from MBD_system.q2dR_i import q2dR_i
from MBD_system.q2theta_i import q2theta_i
from MBD_system.q2dtheta_i import q2dtheta_i
from MBD_system.transform_cs import gcs2cm_lcs
from MBD_system.transform_cs import uP_gcs2lcs
from MBD_system.contact_model.contact_model import ContactModel
from MBD_system.friction_model.friction_model import FrictionModel

class ContactSphereSphere(Contact):
    """
    classdocs
    """
#     __id = itertools.count(0)
    def __init__(self, _type, body_id_i, body_id_j, R0_i, R0_j, properties_dict=[], parent=None):
        #    number
#         self.contact_id = self.__id.next()
        #    name as string
        self._name = "Contact_Sphere_Sphere"# + str(self.contact_id)

        #    this has to be after attributes contact_id and _name as name is constructed from contact_id and _name
        super(ContactSphereSphere, self).__init__(_type, body_id_i, body_id_j, name=self._name, properties_dict=properties_dict, parent=parent)
        """
        Constructor of class contact
        in:
            joint_name - string
            body_i -
            body_j -
        """
        self._parent = parent

        #    contact type
        self._type = "Contact Sphere-Sphere"
        
        #    radius of spheres
        self.R0_i = R0_i
        self.R0_j = R0_j
        self.R0_list = [self.R0_i, self.R0_j]
        
        #    minimum distance value at contact
        self._distance0 = self.R0_i + self.R0_j

        self._contact_point_found = False

        # #   contact model
        # self.contact_model = ContactModel(self.contact_model_type, c_r = self.coef_of_restitution, parent=self)
        #
        # #    friction model
        # self.friction_model = FrictionModel(self.friction_model_type, coef_of_friction_dynamic=self.coef_of_friction_dynamic, coef_of_friction_static=self.coef_of_friction_static, parent=self)

    def check_for_contact(self, q):
        """
        Function check for contact between two spheres
        returns:
                -1 - first contact indentation is too large, return to previous time step and continue with smaller step size
                0 - no contact
                +1 - contact detected
        """
        # print "check_for_contact()"
        #    distance value, distance object is calculated in submethod below
        self._delta = self._evaluate_contact_distance(q)

        #   evaluate normal
        self._n = self._distance_obj.get_normal_2D()

        #   calculate current contact geometry
        self._get_contact_geometry_data(q)

        #   calculate relative contact velocity in normal direction at time t
        self._dq_n, self._dq_t = self._contact_velocity(q)
        # print "self._dq_n =", self._dq_n
        #   add calculated data to container
        if self._step_solution_accepted:
            self._distance_solution_container = np.append(self._distance_solution_container, self._delta)  # _distance_sign, _distance
        else:
            self._distance_solution_container[-1] = self._delta
        
        #   check sign of product between last to values of distance
        self._sign_check = np.sign(self._distance_solution_container[-1]*self._distance_solution_container[-2])
        # print "self._distance =", self._delta
        #    contact
        if self._delta <= 0 and self._dq_n < 0:
            #   condition to find first contact
            if self._sign_check == -1 and self._delta > -self.distance_TOL:
                self._delta0 = self._delta
                self.contact_detected = True
                self.contact_distance_inside_tolerance = True
                self.status = 1
                self._status_container = np.append(self._status_container, self.status)
            #   condition to check if penetration depth of first contact is too large
            else:# self._sign_check == -1 and self._distance > -self.distance_TOL:
                self.status = -1
            # #   no contact
            # else:
            #     self.contact_distance_inside_tolerance = False
            #     self.status = 1
            #     self._status_container = np.append(self._status_container, self.status)
        #    no contact
        if self._delta > 0:
            self.contact_detected = False
            self.contact_distance_inside_tolerance = False
            self.status = 0
            self._status_container = np.append(self._status_container, self.status)
            self.no_contact()
        return self.status

    def contact_update(self, step, t, q):
        """
        Function updates status only for contact type - 2 speheres contact
        :return:
            status - status of contact 0-no contact, 2-contact
        """
        #   step value
        self._step = step

        #   calculate relative contact velocity in normal direction at time t
        self._dq_n, self._dq_t = self._contact_velocity(q)

        #   calculate distance: node-edge
        # self._distance, _inside = evaluate_distance_2D(self.node_GCS, self.edge_GCS[0], self._n, self._t)
        self._delta = self._evaluate_contact_distance(q)

        if self._delta <= self._delta0:#self._delta0:# and self._dq_n <= abs(self._dq0_n):#_inside and self._distance < -1*self.distance_TOL:
            self.status = 1
        else:
            self.status = 0
            self._contact_point_found = False
            self.initial_contact_velocity_calculated = False
            self.contact_distance_inside_tolerance = False
            self.contact_detected = False
            self._status_container = np.append(self._status_container, self.status)

            #   reset values of contact force object to zero if no contact
            _F0 = np.zeros(2)
            _u_P0 = np.zeros(2)
            for force_n, force_t in zip(self._Fn_list, self._Ft_list):
                force_n.update(self._step, F=_F0, u_P=_u_P0)
                force_t.update(self._step, F=_F0, u_P=_u_P0)
        return self.status

    def _get_contact_geometry_data(self, q):
        """
        Function calculates a vector - point of contact from global coordinates to local coordinates of each body
        """
        #   evaluate normal
        self._n = self._distance_obj.get_normal_2D()
        # print "self._n =", self._n
        self._n_list_GCS = [self._n, -self._n]

        #   evaluate tangent
        #   tangent is calculated from rotation of normal for 90deg in CCW direction
        self._t = np.dot(A_matrix(np.pi/2), self._n)
        self._t_list = []

        #   evaluate contact points on each body in body LCS
        self._n_list = []
        self.u_P_list_LCS = []
        for body_id, _normal, _R0 in zip(self.body_id_list, self._n_list_GCS, self.R0_list):
            #   normal in LCS
            _theta = q2theta_i(q, body_id)
            _normal_LCS = uP_gcs2lcs(u_P=_normal, theta=_theta)
            #   append normal to list
            self._n_list.append(_normal_LCS)

            #   tangent in LCS
            _tangent_LCS = np.dot(A_matrix(np.pi/2), _normal)
            self._t_list.append(_tangent_LCS)

            #   calculate actual contact point in LCS on a body surface
            _u_P = _R0 * _normal_LCS

            #   append to contact point u_P vector to list
            self.u_P_list_LCS.append(_u_P)

    def _contact_velocity(self, q):
        """
        Function evaluates relative contact velocity vectors in normal and tangent direction
        """
        dr_P = []
        for body_id, _R, _n, _t in zip(self.body_id_list, self.R0_list, self._n_list, self._t_list):
            #   contact point in body LCS
            _uP = _R * _n

            #   body velocity, R, theta
            _dR = q2dR_i(q, body_id)
            _theta = q2theta_i(q, body_id)

            #    dtheta - omega
            _dtheta = q2dtheta_i(q, body_id)

            #    velocity at contact point (in GCS) on a body
            _dr_P = dr_contact_point_uP(_dR, _theta, _dtheta, _uP)
            
            #    appned to list
            dr_P.append(_dr_P)
        
        #    relative contact velocity vector
        _dq = dr_P[1] - dr_P[0]
        
        #   relative contact velocity (scalar value)
        #   normal direction
        _dq_n = np.dot(_dq, self._n)
        #   tangent direction
        _dq_t = np.dot(_dq, self._t)
        return _dq_n, _dq_t

    def _evaluate_contact_distance(self, q):
        """
        Function evaluates contact distance for this type of contact
        :return:
        """
        #    create distance object
        self._distance_obj = Distance(q2R_i(q, self.body_id_i), q2R_i(q, self.body_id_j))

        #    distance value
        _distance = self._distance_obj._distance - self._distance0
        return _distance

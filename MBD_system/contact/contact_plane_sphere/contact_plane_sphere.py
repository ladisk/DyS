__author__ = 'lskrinjar'

import time
import numpy as np
from pprint import pprint

from MBD_system.contact.contact import Contact
from MBD_system.q2R_i import q2R_i
from MBD_system.contact.distance.distance import Distance
from MBD_system.contact.distance.distance_line_node import DistanceLineNode
from MBD_system.A import A_matrix
from MBD_system.dr_contact_point_uP import dr_contact_point_uP
from MBD_system.q2dR_i import q2dR_i
from MBD_system.q2theta_i import q2theta_i
from MBD_system.q2dtheta_i import q2dtheta_i
from MBD_system.transform_cs import gcs2cm_lcs
from MBD_system.transform_cs import uP_gcs2lcs
from MBD_system.u_P_lcs2gcs import u_P_lcs2gcs
from MBD_system.contact_model.contact_model import ContactModel
from MBD_system.friction_model.friction_model import FrictionModel
from MBD_system.Ai_ui_P import Ai_ui_P_vector
from MBD_system.n2t import n2t


class ContactPlaneSphere(Contact):
    """
    classdocs
    """


    def __init__(self, _type, body_id_i, body_id_j, u_iP, R0_j, n_i, u_jP=np.zeros(2), u_iR=None, t_i=None, properties_dict=[], parent=None):
        """
        Constructor of (sub)class contact of type plane-sphere
        (as we are only working in 2D this is actually a line-circle contact)
        :param body_id_i:       id of line (plane) body
        :param body_id_j:       id of circle (sphere) body
        :param u_iP:            a point vector in LCS of line of body i
        :param R0_j:            radius of a sphere of body j
        :param n_i:             normal vector in LCS of body i
        _param u_iR:            a point vector in LCS of line of body i - optional
        :param t_i:             tangent vector in LCS of i (plane) - optional
        """
        #    name as string
        self._name = "Contact_Plane_Sphere"# + str(self.contact_id)

        #    this has to be after attributes contact_id and _name as name is constructed from contact_id and _name
        super(ContactPlaneSphere, self).__init__(_type, body_id_i, body_id_j, name=self._name, properties_dict=properties_dict, parent=parent)

        self._parent = parent

        #    contact type
        self._type = "Contact Plane-Sphere"
        
        #    radius of contact surfaces
        self.R0_i = None
        self.R0_j = R0_j
        self.R0_list = [self.R0_i, self.R0_j]
        self.u_jP_LCS = u_jP
        self.r_jP_GCS = np.zeros(2)

        #   contact line, plane point
        #   in LCS
        self.u_iP_LCS = u_iP
        self.u_iR_LCS = u_iR
        self.u_i_list_LCS = [self.u_iP_LCS, self.u_iR_LCS]
        #   in GCS
        self.r_i_list_GCS = []

        #   normal of a plane in LCS
        self.n_i_LCS = n_i

        if t_i is None:
            #    calculate tangent based on normal
            self.t_i_LCS = Ai_ui_P_vector(self.n_i_LCS, np.pi/2)
        else:
            self.t_i_LCS = t_i
        
        #    minimum distance value at contact
        self._distance0 = self.R0_j

        #   contact model
        self._create_contact_model(self.properties)

        #   friction model
        self._create_friction_model(self.properties)

        #   list of markers
        self.markers = self._create_markers()

        pprint(vars(self))

    def _create_markers(self):
        """

        :return:
        """
        return []

    def check_for_contact(self, q):
        """
        Function check for contact between plane and sphere
        """
        #   evaluate distance, delta
        self._distance, self._delta = self._contact_geometry_GCS(q)
        # print "---------------------------------"
        # print "self._distance =", self._distance
        # print "self._delta =", self._delta
        #   calculate current contact geometry
        # self._get_contact_geometry_data(q)

        #   calculate relative contact velocity in normal direction at time t
        self._dq_n, self._dq_t = self._contact_velocity(q)

        #   add calculated data to container
        # if self._step_solution_accepted:
        #     self._distance_solution_container.append(self._delta)  # _distance_sign, _distance
        # else:
        #     self._distance_solution_container[-1] = self._delta
        
        #   check sign of product between last to values of distance
        # self._sign_check = np.sign(self._distance_solution_container[-1]*self._distance_solution_container[-2])

        #   check sign
        self._sign_check = np.sign(self._delta * self._distance_solution_container[-1])
        #    contact
        if (self._sign_check == -1) and (self._delta <= 0):
            #   condition to find first contact
            if abs(self._delta) <= self.distance_TOL:
                print "contact detected!"
                # time.sleep(1)
                self._delta0 = self._delta
                self.contact_detected = True
                self.contact_distance_inside_tolerance = True
                self.status = 1

            #   condition to check if penetration depth of first contact is too large
            else:# self._sign_check == -1 and self._distance > -self.distance_TOL:
                self.status = -1
            # #   no contact
            # else:
            #     self.contact_distance_inside_tolerance = False
            #     self.status = 1

        #    no contact
        if self._delta > 0:
            self.contact_detected = False
            self.contact_distance_inside_tolerance = False
            self.status = 0
            self.no_contact()
        return self.status

    def contact_update(self, step, t, q):
        """
        Function updates status only for contact type - 2 speheres contact
        :return:
            status - status of contact 0-no contact, 2-contact
        """
        # print "contact_update()"
        # time.sleep(1)
        #   step value
        self._step = step

        #   calculate relative contact velocity in normal direction at time t
        self._dq_n, self._dq_t = self._contact_velocity(q)

        #   calculate distance: node-edge
        # self._distance, _inside = evaluate_distance_2D(self.node_GCS, self.edge_GCS[0], self._n, self._t)
        self._delta = self._evaluate_contact_distance(q)

        if self._delta <= self._delta0:# and self._dq_n < 0:# and self._dq_n <= abs(self._dq0_n):#_inside and self._distance < -1*self.distance_TOL:
            self.status = 1
        else:
            self.status = 0
            self._contact_point_found = False
            self.initial_contact_velocity_calculated = False
            self.contact_distance_inside_tolerance = False
            self.contact_detected = False
            self.no_contact()


            # #   reset values of contact force object to zero if no contact
            # _F0 = np.zeros(2)
            # _u_P0 = np.zeros(2)
            # for force_n, force_t in zip(self._Fn_list, self._Ft_list):
            #     force_n.update(self._step, F=_F0, u_P=_u_P0)
            #     force_t.update(self._step, F=_F0, u_P=_u_P0)
        return self.status

    def _get_contact_geometry_data(self, q):
        """
        Function calculates a vector - point of contact from global coordinates to local coordinates of each body
        """
        #   evaluate normal for each body in contact
        #   in GCS
        self._n_GCS_list = [+self._n_GCS, -self._n_GCS]

        #   in LCS
        #   list of normals in LCS
        self._n_LCS_list = []
        for body_id, n_i in zip(self.body_id_list, self._n_GCS_list):
            #   normal in LCS
            _theta = q2theta_i(q, body_id)
            _normal_LCS = uP_gcs2lcs(u_P=n_i, theta=_theta)
            #   append normal to list
            self._n_list.append(_normal_LCS)

        # print "self._n_list_GCS =", self._n_list_GCS
        #   evaluate tangent
        #   tangent is calculated from rotation of normal for 90deg in CCW direction
        self._t_GCS = np.dot(A_matrix(np.pi/2), self._n_GCS)
        self._t_GCS_list = [self._t_GCS, -self._t_GCS]

        #   contact point on body j in GCS
        self.r_jP_GCS = q2R_i(q, self.body_id_j) + self.R0_j * self._n_GCS_list[1]

        #   contact point on body i in GCS
        u_iP_LCS = self._distance_obj.contact_point_on_line() + self.u_iP_LCS
        self.r_iP_GCS = u_P_lcs2gcs(u_iP_LCS, q, self. body_id_i)

        #   list of contact point in GCS
        self.u_P_GCS_list = [self.r_iP_GCS, self.r_jP_GCS]


    # def _contact_velocity(self, q):
    #     """
    #     Function evaluates relative contact velocity vectors in normal and tangent direction
    #     """
    #     dr_P = []
    #     for body_id, _R, _n, _t in zip(self.body_id_list, self.R_list, self._n_list, self._t_list):
    #         #   contact point in body LCS
    #         _uP = _R * _n
    #
    #         #   body velocity, R, theta
    #         _dR = q2dR_i(q, body_id)
    #         _theta = q2theta_i(q, body_id)
    #
    #         #    dtheta - omega
    #         _dtheta = q2dtheta_i(q, body_id)
    #
    #         #    velocity at contact point (in GCS) on a body
    #         _dr_P = dr_contact_point_uP(_dR, _theta, _dtheta, _uP)
    #
    #         #    appned to list
    #         dr_P.append(_dr_P)
    #
    #     #    relative contact velocity vector
    #     _dq = dr_P[1] - dr_P[0]
    #
    #     #   relative contact velocity (scalar value)
    #     #   normal direction
    #     _dq_n = np.dot(_dq, self._n)
    #     #   tangent direction
    #     _dq_t = np.dot(_dq, self._t)
    #     return _dq_n, _dq_t

    def _evaluate_contact_distance(self, q):
        """
        Function evaluates contact distance for this type of contact
        :return:
        """
        self._contact_geometry_GCS(q)
        #    create distance object
        # print "q2R_i(q, self.body_id_i) =", q2R_i(q, self.body_id_i)
        # print "q2R_i(q, self.body_id_j) =", q2R_i(q, self.body_id_j)
        #   q2R_i(q, self.body_id_i), q2R_i(q, self.body_id_j)
        # self._distance_obj = Distance(q2R_i(q, self.body_id_j), self.u_i_list[0], n2=self.u_i_list[1], normal=self.n_i)
        self._distance_obj = DistanceLineNode(q2R_i(q, self.body_id_j), self.r_P_GCS_list[0], self._n_GCS)
        # pprint(vars(self._distance_obj))
        #    distance value
        _distance = self._distance_obj._distance - self._distance0

        # print "_distance(_evaluate_contact_distance) =", _distance
        # time.sleep(10)
        return _distance

    def _contact_geometry_GCS(self, q):
        """
        Function evaluates contact geometry data in GCS
        :param q:
        :return:
        """
        self.r_i_list_GCS, self.r_jP_GCS = self._contact_geometry_P_GCS(q)

        #   plane normal in GCS
        self._n_GCS = Ai_ui_P_vector(self.n_i_LCS, q2theta_i(q, self.body_id_i))

        #   distance object
        self._distance_obj = DistanceLineNode(q2R_i(q, self.body_id_j), self.r_P_GCS_list[0], normal=self._n_GCS)

        #   get distance and delta value
        _distance = self._distance_obj._distance

        _delta = _distance - self.R0_j

        return _distance, _delta

    def _contact_geometry_P_GCS(self, q):
        """
        Function evaluates position of body points (sphere center) and points that define edge - line
        :return:
        """
        r_i_list_GCS = []
        for i, body_id in enumerate(self.body_id_list):
            #   geometry points of line
            if i == 0:
                for u_iP in self.u_i_list_LCS:
                    if u_iP is not None:
                        r_iP = u_P_lcs2gcs(u_iP, q, body_id)
                        r_i_list_GCS.append(r_iP)
            #   circle center
            if i == 1:
                r_jP_GCS = u_P_lcs2gcs(self.u_jP_LCS, q, body_id)

        return r_i_list_GCS, r_jP_GCS

    def _contact_geometry_LCS(self, q):
        """
        Function calculates the contact geometry from GCS to LCS (only once).
        Especially uPi, uPj
        :param q:
        :return:
        """
        #   predefine (empty) list of normals in LCS
        self._n_LCS_list = []

        #   evaluate contact points on each body in body LCS
        self.u_P_LCS_list = []
        for body_id, r_P, _normal in zip(self.body_id_list, self.u_P_GCS_list, self._n_GCS_list):
            # uP_lcs = self.u_P_GCS - q2R_i(q, body_id)

            # _theta = q2theta_i(q, body_id)
            # _u_P = uP_gcs2lcs(u_P=uP_lcs, theta=_theta)
            u_P = gcs2cm_lcs(r_P, q2R_i(q, body_id), q2theta_i(q, body_id))
            #   append to contact point u_P vector to list
            self.u_P_LCS_list.append(u_P)

    def _paint_GL_GCS(self, step=None):
        """
        Paint contact point in GCS
        :return:
        """
        #   during integration process (simulation)
        if step is None:
            #   display contact point in GCS on each body in contact
            if self.contact_detected:
                pass

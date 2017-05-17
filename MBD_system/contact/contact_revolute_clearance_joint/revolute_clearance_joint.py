"""
Created on 21. feb. 2014

@author: lskrinjar
"""
import time
import copy
from pprint import pprint
import numpy as np
from OpenGL.GL import *


from MBD_system.A import A_matrix
from MBD_system.contact.contact import Contact
from MBD_system.contact.distance.distance_RCJ import DistanceRCJ
from MBD_system.contact_model.contact_model_cylinder import ContactModelCylinder
from MBD_system.dr_contact_point_uP import dr_contact_point_uP
from MBD_system.q2dR_i import q2dR_i
from MBD_system.q2dtheta_i import q2dtheta_i
from MBD_system.q2q_body import q2q_body
from MBD_system.q2theta_i import q2theta_i
from MBD_system.transform_cs import gcs2cm_lcs
from MBD_system.transform_cs import uP_gcs2lcs
from MBD_system.u_P_lcs2gcs import u_P_lcs2gcs
from MBD_system.q2R_i import q2R_i
from MBD_system.n2t import n2t


class RevoluteClearanceJoint(Contact):
    """
    classdocs
    """


    def __init__(self, _type, body_id_i, body_id_j, u_iP, u_jP, R0_i, R0_j, properties_dict=[], parent=None):
        """
        Constructor of class contact of revolute clearance joint
        :param _type:       type of clearance joint or contact
        :param body_id_i:   id of body - hole (bearing)
        :param body_id_j:   id of body - pin (journal) 
        :param u_iP:        vector to center of a hole on body i in body LCS
        :param u_jP:        vector to center of a pin on body j in body LCS
        :param R0_i:        radius of a hole
        :param R0_j:        radius of a pin
        :param properties_dict: additioanl parameters to override default values or add new parameters
        """
        #    name as string
        self._name = "RC_Joint_"# + str(self.contact_id)

        #    this has to be after attributes contact_id and _name as name is constructed from contact_id and _name
        super(RevoluteClearanceJoint, self).__init__(_type, body_id_i, body_id_j, name=self._name, properties_dict=properties_dict, parent=parent)

        self._parent = parent

        #    type of contact
        self._contact_type = "revolute clearance joint"

        #   vector of axis on revolute joint in LCS of a body i, j
        self.u_iP = u_iP
        self.u_jP = u_jP
        #   predefined empty list of center point or clearance joint (to center of pin/hole) in body LCS
        self.u_CP_LCS_list = [self.u_iP, self.u_jP]
        #   centers of revolute clearance joint in GCS
        self.u_iCP_GCS = np.zeros(2)
        self.u_jCP_GCS = np.zeros(2)
        self.r_CP_GCS_list = [self.u_iCP_GCS, self.u_jCP_GCS]

        #   clearance parameters
        self.R0_i = R0_i
        self.R0_j = R0_j
        self.R0_list = [self.R0_i, self.R0_j]
        #   radial clearance
        self._radial_clearance = abs(R0_i - R0_j)

#         if self.contact_model_type is None:
#             self.contact_model_type = "hertz"
#         # print "self.properties_contact_model =", self.properties_contact_model
#         self.contact_model = ContactModelCylinder(self.contact_model_type, properties_dict=self.properties_contact_model, parent=self)

        #    joint body ids
        self.body_id_i = body_id_i
        self.body_id_j = body_id_j
        self.body_id_list = [self.body_id_i, self.body_id_j]
        
        #   contact model
        self._create_contact_model(self.properties)
        
        #   friction model
        self._create_friction_model(self.properties)

        #   list of markers
        self.markers = self._create_markers()

    def _solution_containers_additional(self):
        """
        Function creates additional solution containers that are specific for type of contact
        :return:
        """
        #   eccentricity vector
        self._e_solution_container = [0, 0]

    def _create_markers(self):
        """
        Function creates markers
        :return:
        """
        markers = []
        for body, _u_P in zip(self.body_list, self.u_CP_LCS_list):
            #    node coordinates
            node = np.array(np.append(_u_P, self.z_dim), dtype='float32')
            
            #    create marker object
            _marker = Marker(node, visible=True, parent=body)
            
            #    append marker to list of body markers
            body.markers.append(_marker)
            #    append marker to list of revolute clearance joint markers
            markers.append(_marker)

        return markers
    
    def check_for_contact_started_condition(self, delta, sign):
        """
        Function checks condition for contact
        """
        if (sign == -1) and (delta <= 0):
            return True
        elif (sign == +1) and (delta <= 0) and (np.sign(self._distance_solution_container[-1]) == -1) and (self._distance_solution_container[-1] > delta):
            return True
        else:
            return False

    def check_for_contact_continued_condition(self, delta, dq_n):
        """
        Function checks condition for contact to continue
        :param delta:
        :return:
        """
        return (delta <= self._delta0)# and (dq_n < self._dq0_n)

    def check_for_contact(self, step, t, q):
        """
        Function checks for contact of revolute clearance joint
        """
        #   step
        self._step = step
        #   time
        self._t = t
        
        #   evaluate distance, delta
        self._distance, self._delta, self._n_GCS, self._t_GCS = self._contact_geometry_GCS(q)

        # self._get_contact_geometry_data(q)
        #   add distance value to container
        # self._distance_solution_container = np.append(self._distance_solution_container, self._delta)

        #   check sign
        self._sign_check = np.sign(self._delta * self._distance_solution_container[-1])

        #    contact has happened, but time step has to be reduced as initial penetration depth is too large
        # if (np.sign(self._dqn_solution_container[-1]) == +1) or (self._dqn_solution_container[-1] == 0) and (self._sign_check == -1) and (self._distance >= self._radial_clearance):
        # if (self._sign_check == -1s) and (self._distance >= self._radial_clearance):
        # if (self._delta <= 0):
        if self.check_for_contact_started_condition(self._delta, self._sign_check):
        # if (self._sign_check == -1) and (self._delta <= 0):
            #    beginning of contact detected, all parameters are within limits
            if abs(self._delta) <= self.distance_TOL:
                # print "check 2"
                self._delta0 = self._delta
                # print "------------------------------"
                # print "contact detected"
                # print "step =", self._step
                # print "self._sign_check =", self._sign_check
                # print "self._distance_solution_container[-1] =", self._distance_solution_container[-1]
                # print "self._delta0 =", self._delta0
                # print "self._n_GCS =", self._n_GCS
                # print "self._t_GCS =", self._t_GCS
                self.contact_detected = True
                self.contact_distance_inside_tolerance = True
                # print "self.contact_detected =", self.contact_detected
                self.status = 1
                return 1

            #   step back
            if abs(self._delta) > self.distance_TOL:
                self.contact_detected = True
                self.status = -1
                # self._delta = 1E-100
                return -1

        #    all calculated distances are greater than tolerance and bodies are not in contact
        self.status = 0
        self.no_contact()

        return 0

    def _contact_update_condition(self):
        """
        Function checks condition for contact if contact is still present
        """
        return (self._delta <= self._delta0) and (self._dq_n < abs(self._dq0_n))

    def contact_update(self, step, t, q):
        """
        Function evaluates contact distance while, contact is present
        :return:
            status - status of contact 0-no contact, 2-contact
        """
        # print "contact_update()"
        self._step = step
        # print "self._step =", self._step
        #   current contact velocity at time t
        self._dq_n, self._dq_t = self._contact_velocity(q)
#         print "self._dq_n =", self._dq_n
        #   calculate distance between joint centers and delta (penetration)
        self._distance, self._delta, self._n_GCS, self._t_GCS = self._contact_geometry_GCS(q)
        # print "self._distance", self._distance
        # print "self._delta", self._delta
        # print "self._distance >= self._radial_clearance =", self._distance >= self._radial_clearance
        # print "abs(self._delta) >= self.distance_TOL =", abs(self._delta) >= self.distance_TOL
        # time.sleep(1)
        #   if distance is greater than radial clearance, contact is present
        if (self._delta <= self._delta0):# and (self._dq_n < abs(self._dq0_n)):#(self._distance >= self._radial_clearance) and (abs(self._delta) >= self.distance_TOL):
            self.status = 1

        else:
            # print "contact finished"

            self.status = 0
            self.contact_detected = False
            self.contact_distance_inside_tolerance = False
            # self.list_of_contact_force_objects_constructed = False
            self._contact_point_found = False
            self.initial_contact_velocity_calculated = False
#             print "self._dq_n =", self._dq_n
            self.no_contact()
#             print "no contact: t =", t, "step =", step
#             print "self._delta0 =", self._delta0
#             print "self._delta =", self._delta

        # print "self.status =", self.status
        return self.status

    def _get_contact_geometry_data(self, q):
        """
        Function calculates normal and tangent vector of body i, j in GCS
        """
        #   normal in GCS of body i, j
        self._n_GCS_list = [+self._n_GCS, -self._n_GCS]

        #   tangent in GCS of body i, j
        self._t_GCS_list = [+self._t_GCS, -self._t_GCS]

    def _contact_geometry_CP_GCS(self, q):
        """
        Function calculates position of centers (CP - Center Points) of revolute joint pin/hole in GCS
        :param q:               a vector of coordinates of the system (numpy array)
        :return u_CP_list_GCS:  a list of two arrays (vectors) of
        """
        #   calculate position of pin/hole centers of each body in GCS
        for i, (body_id, u_P) in enumerate(zip(self.body_id_list, self.u_CP_LCS_list)):
            #   axis center of revolute joint of each body in GCS
            self.r_CP_GCS_list[i] = u_P_lcs2gcs(u_P, q, body_id)

        #   convert to 32bit float to display in opengl scene
        return np.array(self.r_CP_GCS_list, dtype="float32")

    def _contact_geometry_GCS(self, q):
        """
        Function calculates distance between center points (eccentricity vector) and indentation based on radius of pin/hole
        Function calculates actual contact point on each body in body LCS and GCS
        :param q:               vector of states of MBD system
        :return _distance:      distance between center of pin and hole
        :paramm _delta:         depth of deformation (indentation, penetration)
        :return _n_GCS:         normal of contact in GCS
        _return _t_GCS:         tangent of contact in GCS
        """
        #   a list of center point of pin and hole on each body in GCS
        self.r_CP_GCS_list = self._contact_geometry_CP_GCS(q)

        #   calculate distance between axis of both bodies in revolute joint
        self._distance_obj = DistanceRCJ(self.r_CP_GCS_list[0], self.r_CP_GCS_list[1], parent=self)

        #   penetration depth is difference between nominal radial clearance and actual calculated clearance at time t
        _distance = self._distance_obj._distance
        _delta = self._radial_clearance - _distance
        
        #   unit vector in direction from one center to another (pin to hole)
        _n_GCS = self._distance_obj._distance_vector / self._distance_obj._distance

        #    tangent in GCS
        _t_GCS = n2t(_n_GCS)

        #   create normal list in GCS
#         self._n_GCS_list = [-self._n_GCS, +self._n_GCS]
#         self._t_GCS_list = []

#         #   calculate a actual contact point in revolute clearance joint of each body in GCS
        self.r_P_GCS_list = copy.copy(self.r_P_GCS_0_list)
        #   evaluate actual contact point in LCS of each body and in GCS
        for i, (_u_CP_GCS, _u_CP_LCS, _R0) in enumerate(zip(self.r_CP_GCS_list, self.u_CP_LCS_list, self.R0_list)):
            #   contact point in GCS
            _u_P_GCS = _u_CP_GCS + _R0 * self._n_GCS
 
            #   add to list
            self.r_P_GCS_list[i] = _u_P_GCS
# 
#             #   create normal in GCS of a body
#             t_i = n2t(n_i)
#             self._t_GCS_list.append(t_i)

        return _distance, _delta, _n_GCS, _t_GCS
    
    def _contact_geometry_LCS(self, q):
        """
        Function evaluates contact geometry parameters in body LCS based on contact geometry in GCS
        Function evaluates:
        self._n_LCS_list:   normal of contact in body LCS:
        self._t_LCS_list:   tangent of contact in body LCS
        self.u_P_LCS_list:  contact point in body LCS
        """
        #   predefine (empty) list of normals in LCS
        self.u_P_LCS_list = []
        self._n_LCS_list = []

        #    vector of contact point in LCS
        for i, (body_id, _normal, r_P) in enumerate(zip(self.body_id_list, self._n_GCS_list, self.r_P_GCS_list)):
            #    R
            R_i = q2R_i(q, body_id)
            #   theta
            theta_i = q2theta_i(q, body_id)

            #   normal in LCS
            normal_LCS = uP_gcs2lcs(u_P=_normal, theta=theta_i)
            #   append normal to list
            self._n_LCS_list.append(normal_LCS)

            #   contact point in body LCS
            u_P = gcs2cm_lcs(r_P, R_i, theta_i)
            #   append to list
            self.u_P_LCS_list.append(u_P)

    def _contact_velocity(self, q):
        """
        Function evaluates relative contact velocity at contact point in revolute clearance joint.
        :param q:   array of coordinates (r, theta) and velocities (dR, dheta=omega) of MBD system
        """
        dr_P = []
        for body_id, u_P in zip(self.body_id_list, self.u_P_LCS_list):
            dR = q2dR_i(q, body_id)
            theta = q2theta_i(q, body_id)
            #    dtheta - omega
            dtheta = q2dtheta_i(q, body_id)

            #    point velocity
            dr_P_body = dr_contact_point_uP(dR, theta, dtheta, u_P)

            #    add to list
            dr_P.append(dr_P_body)

        _dq = dr_P[1] - dr_P[0]

        #   relative contact velocity
        #   normal direction
        _dq_n = np.dot(_dq, self._n_GCS)

        #   tangent direction
        _dq_t = np.dot(_dq, self._t_GCS)

        return _dq_n, _dq_t

#     def solve(self, t, q):
#         """
#         Calculate contact parameters
#         returns:
#         """
#         # print "solve()"
#         #    calculate coordinates of contact point from global coordinates in local coordinates of each body in contact
#         if not self._contact_point_found:
#             self._get_contact_geometry_data(q)
#             self._contact_point_found = True
#         else:
#             pass
#             # self._distance, self._delta = self._contact_geometry_GCS(q)
#             # print "contact update"
# 
#         # print "self._contact_point_found =", self._contact_point_found
# 
#         #   kinematic properties of contact point
#         #   initial contact velocity
#         if not self.initial_contact_velocity_calculated:
#             self._dq0_n, self._dq0_t = self._contact_velocity(q)
#             self.contact_model.set_dq0(self._dq0_n, self._dq0_t)
#             self.initial_contact_velocity_calculated = True
# 
#         # if self.contact_detected:
#         self._distance, self._delta = self._contact_geometry_GCS(q)
#         # print "self._delta =", self._delta
#         # print "self._dq0_n, self._dq0_t =", self._dq0_n, self._dq0_t
# 
#         #   current contact velocity at time t
#         self._dq_n, self._dq_t = self._contact_velocity(q)
#         # print "self._dq_n, self._dq_t =", self._dq_n, self._dq_t
#         # time.sleep(100)
#         #   current contact velocity at time t
#         # self._dq_n, self._dq_t = self._contact_velocity(q)
# 
#         if self._contact_type.lower() in self._contact_types:#== "general" or self._type.lower() == "revolute clearance joint" or self._type.lower() == "contact sphere-sphere":#ECF-N
#             self._solve_ECF_N(t, q, self._delta, self._dq_n, self._dq_t)#self._delta
#         else:
#             raise ValueError, "Contact type not correct!"
        
    def evaluate_rijP(self, q):
        """
        Evaluate distance vector - eccentricity vector
        """
        self.r_CP_GCS_list = self._contact_geometry_CP_GCS(q)
        #   calculate distance between axis of both bodies in revolute joint
        self._distance_obj = DistanceRCJ(self.r_CP_GCS_list[0], self.r_CP_GCS_list[1], parent=self)

        #   penetration depth is difference between nominal radial clearance and actual calculated clearance at time t
        _distance = self._distance_obj._distance
        
        return _distance
    
    def evaluate_delta(self, q):
        """
        
        """
        self.r_CP_GCS_list = self._contact_geometry_CP_GCS(q)
        #   calculate distance between axis of both bodies in revolute joint
        self._distance_obj = DistanceRCJ(self.r_CP_GCS_list[0], self.r_CP_GCS_list[1], parent=self)

        #   penetration depth is difference between nominal radial clearance and actual calculated clearance at time t
        _distance = self._distance_obj._distance
        
        #    delta
        _delta = self._radial_clearance - _distance
        
        return _delta

    def _paint_GL_GCS(self, step=None):
        """
        Paint contact point in GCS
        :return:
        """
        #   during integration process (simulation)
        if step is None:
            #   display contact point in GCS on each body in contact
            if self.contact_detected:

                #   actual contact point on body i, j in GCS
                for body_id, r_P_GCS, Fn, Ft in zip(self.body_id_list, self.r_P_GCS_list, self._Fn_list, self._Ft_list):
                    if Fn._visible or Ft._visible:
                        glColor3f(self._bodies[body_id].color_GL[0], self._bodies[body_id].color_GL[1], self._bodies[body_id].color_GL[2])
                        glVertex3f(r_P_GCS[0], r_P_GCS[1], self.z_dim)

            #   center of pin and hole
            if self._distance_obj is not None:
                for body_id, r_P in zip(self.body_id_list, [self._distance_obj.r_iP, self._distance_obj.r_jP]):
                    #    paint pin center in GCS
                    glColor3f(self._bodies[body_id].color_GL[0], self._bodies[body_id].color_GL[1], self._bodies[body_id].color_GL[2])
                    glVertex3f(r_P[0], r_P[1], self.z_dim)

        #   during animation
        else:
            #   paint contact point only if contact is present (status=1)
            if self._status_container[step] == 1:
                for i, (body_id, Fn, Ft) in enumerate(zip(self.body_id_list, self._Fn_list, self._Ft_list)):
                    if Fn._visible or Ft._visible:
                        glColor3f(self._bodies[body_id].color_GL[0], self._bodies[body_id].color_GL[1], self._bodies[body_id].color_GL[2])
                        glVertex3f(self._r_P_solution_container[step][i][0], self._r_P_solution_container[step][i][1], self.z_dim)

    def testing(self):
#         #    evaluate rCP vector of center points in GCS of pin and hole
#         q = self._parent._parent.evaluate_q()
#         print self._contact_geometry_CP_GCS(q)
        print self.contact_model
        pprint(vars(self.contact_model))


if __name__ == "__main__":
    a = RevoluteClearanceJoint(body_id_i=0, body_id_j=1, u_iP=np.array([0, 0]), u_jP=np.array([-1, 0]), type="rcj", properties_dict=[], parent=None)
    print "a =", a
    pprint(vars(a))
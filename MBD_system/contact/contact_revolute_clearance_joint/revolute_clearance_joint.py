"""
Created on 21. feb. 2014

@author: lskrinjar
"""
from pprint import pprint

import numpy as np

from MBD_system.A import A_matrix
from MBD_system.contact.contact import Contact
from MBD_system.contact.distance.distance_revolute_clearance_joint import DistanceRevoluteClearanceJoint
from MBD_system.contact_model.contact_model_cylinder import ContactModelCylinder
from MBD_system.dr_contact_point_uP import dr_contact_point_uP
from MBD_system.q2dR_i import q2dR_i
from MBD_system.q2dtheta_i import q2dtheta_i
from MBD_system.q2q_body import q2q_body
from MBD_system.q2theta_i import q2theta_i
from MBD_system.transform_cs import gcs2cm_lcs
from MBD_system.transform_cs import uP_gcs2lcs
from MBD_system.u_P_lcs2gcs import u_P_lcs2gcs
from simulation_control_widget.opengl_widget.marker.marker import Marker


class RevoluteClearanceJoint(Contact):
    """
    classdocs
    """
#     __id = itertools.count(0)
    def __init__(self, _type, body_id_i, body_id_j, u_iP, u_jP, R0_i, R0_j, properties_dict=[], parent=None):
        """
        Constructor of class contact of revolute clearance joint
        :param _type:       type of clearance joint or contact
        :param body_id_i:   id of hole body
        :param body_id_j:   id of pin body
        :param u_iP:        vector to center of a hole on body i in body LCS
        :param u_jP:        vector to center of a pin on body j in body LCS
        :param R0_i:        radius of a hole
        :param R0_j:        radius of a pin
        :param properties_dict: additioanl parameters to override default values or add new parameters
        """
        #    number
#         self.contact_id = self.__id.next()
        #    name as string
        self._name = "RC_Joint_"# + str(self.contact_id)

        #    this has to be after attributes contact_id and _name as name is constructed from contact_id and _name
        super(RevoluteClearanceJoint, self).__init__(_type, body_id_i, body_id_j, name=self._name, properties_dict=properties_dict, parent=parent)

        self._parent = parent

        #    joint type
        self._type = "Revolute Clearance Joint"

        #   vector of axis on revolute joint in LCS of a body i, j
        self.u_iP = u_iP
        self.u_jP = u_jP
        #   predefined empty list of center point or clearance joint (to center of pin/hole) in body LCS
        self.u_CP_LCS_list = [self.u_iP, self.u_jP]
        #   centers of revolute clearance joint in GCS
        self.u_iCP_GCS = np.zeros(2)
        self.u_jCP_GCS = np.zeros(2)
        self.u_CP_GCS_list = [self.u_iCP_GCS, self.u_jCP_GCS]

        #   clearance parameters
        self.R0_i = R0_i
        self.R0_j = R0_j
        self.R0_list = [self.R0_i, self.R0_j]
        #   radial clearance
        self._radial_clearance = abs(R0_i - R0_j)

        #   contact model
        if self.contact_model_type is None:
            self.contact_model_type = "hertz"
        # print "self.properties_contact_model =", self.properties_contact_model
        self.contact_model = ContactModelCylinder(self.contact_model_type, properties_dict=self.properties_contact_model, parent=self)

        #    joint body ids
        self.body_id_i = body_id_i
        self.body_id_j = body_id_j
        self.body_id_list = [self.body_id_i, self.body_id_j]

        #   list of markers
        self.markers = self._create_markers()

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

    def check_for_contact(self, q):
        """
        Function check for contact of revolute clearance joint
        """
        # print "check_for_contact()"
        #   evaluated distance, delta
        self._distance, self._delta = self._contact_geometry_GCS(q)
        print "self._delta =", self._delta
        #   add distance value to container

        # self._distance_solution_container = np.append(self._distance_solution_container, self._delta)

        #   check sign
        self._sign_check = np.sign(self._delta * self._distance_solution_container[-1])

        #    contact has happened, but time step has to be reduced as initial penetration depth is too large
        # if (np.sign(self._dqn_solution_container[-1]) == +1) or (self._dqn_solution_container[-1] == 0) and (self._sign_check == -1) and (self._distance >= self._radial_clearance):
        # if (self._sign_check == -1) and (self._distance >= self._radial_clearance):
        if (self._sign_check == -1) and (self._delta <= 0):

            #    beginning of contact detected, all parameters are within limits
            if abs(self._delta) <= self.distance_TOL:
            # else:
                self._delta0 = self._delta
                self.contact_detected = True
                self.contact_distance_inside_tolerance = True
                self.status = 1
                # print int(self._step_num_solution_container[-1]), self.status, t, self._distance,
                return 1

            #   step back
            if abs(self._delta) > self.distance_TOL:
                self.contact_detected = True
                self.status = -1
                # self._distance_solution_container = np.delete(self._distance_solution_container, -1)
                return -1

        #    all calculated distances are greater than tolerance and bodies are not in contact
        self.status = 0
        # self._status_container = np.append(self._status_container, self.status)
        self.no_contact()

        # for _force_n, _force_t in zip(self._Fn_list, self._Ft_list):
        #     _force_n.update(self._step)
        #     _force_t.update(self._step)

        return 0

    def contact_update(self, step, t, q):
        """
        Function evaluates contact distance while, contact is present
        :return:
            status - status of contact 0-no contact, 2-contact
        """
        self._step = step
        # print "self._step =", self._step
        #   current contact velocity at time t
        self._dq_n, self._dq_t = self._contact_velocity(q)
        # print "self._dq_n =", self._dq_n
        #   calculate distance between joint centers and delta (penetration)
        self._distance, self._delta = self._contact_geometry_GCS(q)
        # print "self._distance", self._distance
        # print "self._delta", self._delta
        # print "self._distance >= self._radial_clearance =", self._distance >= self._radial_clearance
        # print "abs(self._delta) >= self.distance_TOL =", abs(self._delta) >= self.distance_TOL
        # time.sleep(1)
        #   if distance is greater than radial clearance, contact is present
        if self._delta <= self._delta0:#(self._distance >= self._radial_clearance) and (abs(self._delta) >= self.distance_TOL):
            self.status = 1

        else:
            self.status = 0
            self.contact_detected = False
            self.contact_distance_inside_tolerance = False
            # self.list_of_contact_force_objects_constructed = False
            self._contact_point_found = False
            self.initial_contact_velocity_calculated = False

            self.no_contact()

        return self.status

    def _get_contact_geometry_data(self, q):
        """
        Function calculates a vector - point of contact from global coordinates to local coordinates of each body
        """
        #   tangent is calculated from rotation of normal for 90deg in CCW direction
        self._t_GCS = np.dot(A_matrix(np.pi/2), self._n_GCS)
        self._t_LCS_list = [self._t_GCS, -self._t_GCS]

    def _contact_geometry_CP_GCS(self, q):
        """
        Function calculates position of centers (CP - Center Points) of revolute joint pin/hole in GCS
        :param q:               a vector of coordinates of the system (numpy array)
        :return u_CP_list_GCS:  a list of two arrays (vectors) of
        """
        #   calculate position of pin/hole centers of each body in GCS
        u_CP_GCS_list = []
        print "q =", q
        for body_id, u_P in zip(self.body_id_list, self.u_CP_LCS_list):
            print "u_P =", u_P
            #   axis center of revolute joint of each body in GCS
            u_P_GCS = u_P_lcs2gcs(u_P, q, body_id)
            print "u_P_GCS =", u_P_GCS
            u_CP_GCS_list.append(u_P_GCS)

        #   reformat to 32bit float to display in opengl scene
        return np.array(u_CP_GCS_list, dtype="float32")

    def _contact_geometry_GCS(self, q):
        """
        Function calculates distance between centerpoints and indentation based on radius of pin/hole
        :param q:
        :return distance_obj:   distance object of calculated data in GCS
        """
        self.u_CP_GCS_list = self._contact_geometry_CP_GCS(q)
        print "self.u_CP_GCS_list =", self.u_CP_GCS_list
        #   calculate distance between axis of both bodies in revolute joint
        self._distance_obj = DistanceRevoluteClearanceJoint(self.u_CP_GCS_list[0], self.u_CP_GCS_list[1], parent=self)

        #   penetration depth is difference between nominal radial clearance and actual calculated clearance at time t
        _distance = self._distance_obj._distance
        _delta = self._radial_clearance - _distance

        #   contact is present and actual contact point has to be evaluated
        # if _distance >= self._radial_clearance and abs(_delta) >= self.distance_TOL:
        # print "contact is present"
        #   unit vector in direction from one center to enother (pin to hole)
        self._n_GCS = self._distance_obj._distance_vector / self._distance_obj._distance
        # print "--------------------------------"
        # print "step =", self._step
        # print "n =", self._n
        # if _delta < 0 and abs(_delta) > self.distance_TOL:
        #     print "body i =", self.u_CP_list_GCS[0]
        #     print "body j =", self.u_CP_list_GCS[1]
        #     print "n =", self._n

        #   create normal list in LCS
        self._n_GCS_list = [-self._n_GCS, +self._n_GCS]
        self._n_LCS_list = []
        for body_id, _normal in zip(self.body_id_list, self._n_GCS_list):
            #   normal in LCS
            _theta = q2theta_i(q, body_id)
            _normal_LCS = uP_gcs2lcs(u_P=_normal, theta=_theta)
            #   append normal to list
            self._n_LCS_list.append(_normal_LCS)

        #   calculate a actual contact point in revolute clearance joint of each body in GCS
        #   and vector of contact point in LCS
        self.u_P_GCS_list = []
        self.u_P_LCS_list = []
        # plot = False#[False, self.status==2] #False
        # if plot:
        #     print "*************************"
        #     fig = plt.gcf()
        #     plt.xlim([0.0, 0.01])
        #     plt.ylim([0.0, 0.01])
        #     ax = plt.gca()
        #     # ax.relim()
        #     ax.autoscale_view(True,True,True)
        #     ax.set_aspect('equal')

        self.u_P_LCS_list = []
        #   evaluate actual contact point in LCS of each body and in GCS
        for body_id, _u_CP_GCS, _u_CP_LCS, _R0 in zip(self.body_id_list, self.u_CP_GCS_list, self.u_CP_LCS_list, self.R0_list):
            # print "body_id =", body_id
            #   q of a body
            q_body = q2q_body(q, body_id)

            #   contact point in GCS
            _u_P_GCS = _u_CP_GCS + _R0 * self._n_GCS

            #   contact point in body LCS
            _u_P_LCS = gcs2cm_lcs(_u_P_GCS, q_body[0:2], q_body[2])
            self.u_P_LCS_list.append(_u_P_LCS)

            # if plot:
            #     print "----------------------------------"
            #     print "body =", self._parent._parent.bodies[body_id]._name
            #
            #     R = q_body[0:2]
            #     print "q_body[0:2] =", q_body[0:2]
            #     print "joint center =", _u_CP_LCS
            #     print "contact point GCS =", _u_P_GCS
            #     print "contact point LCS =", _u_P_LCS
            #     _color = np.random.rand(3)
            #     circle=plt.Circle((_u_CP_LCS[0]+R[0],_u_CP_LCS[1]+R[1]),_R,color=_color, fill=False)
            #     #   center of axis
            #     ax.plot(_u_CP_LCS[0], _u_CP_LCS[1], "o", color=_color)
            #     #   contact point
            #     ax.plot(_u_P_LCS[0]+R[0], _u_P_LCS[1]+R[1], "x", color=_color)
            #     #   LCS
            #     ax.plot(R[0], R[1], "s", color=_color)
            #     fig.gca().add_artist(circle)

        #   transform to 32bit float to display in opengl
        self.u_P_GCS_list = np.array(self.u_P_GCS_list, dtype="float32")
        self.u_P_LCS_list = np.array(self.u_P_LCS_list, dtype="float32")
        # self.u_P_GCS = np.array(_u_P_GCS, dtype="float32")

        # if self._step_solution_accepted:
        #     self._u_P_solution_container.append(self.u_P_list_LCS.flatten())

        # if plot:
        #     plt.grid()
        #     plt.show()
        #     fig.savefig('testing.png')

        return _distance, _delta

    def _contact_velocity(self, q):
        """
        Function evaluates relative contact velocity at contact point in revolute clearance joint.
        :param q:   array of coordinates (r, theta) and velocities (dR, dhete=omega) of MBD system
        """
        dr_P = []
        print "self.u_P_LCS_list =", self.u_P_LCS_list, type(self.u_P_LCS_list)
        for body_id, u_P in zip(self.body_id_list, self.u_P_LCS_list):
            dR = q2dR_i(q, body_id)
            theta = q2theta_i(q, body_id)
            #    dtheta - omega
            dtheta = q2dtheta_i(q, body_id)

            #    point velocity
            dr_P_body = dr_contact_point_uP(dR, theta, dtheta, u_P)

            #    add to list
            dr_P.append(dr_P_body)

        _dq = dr_P[0] - dr_P[1]

        #   relative contact velocity
        #   normal direction
        _dq_n = np.dot(_dq, self._n_GCS)

        #   tangent direction
        _dq_t = np.dot(_dq, self._t_GCS)

        return _dq_n, _dq_t

    def solve(self, t, q):
        """
        Calculate contact parameters
        returns:
        """
        # print "solve() @ RCJ()"
        # print "Ry_i =", q[1]
        # print "Ry_j =", q[4]
        #    calculate coordinates of contact point from global coordinates in local coordinates of each body in contact
        if not self._contact_point_found:
            self._get_contact_geometry_data(q)
            self._contact_point_found = True
        else:
            pass
            # self._distance, self._delta = self._contact_geometry_GCS(q)
            # print "contact update"

        # print "self._contact_point_found =", self._contact_point_found

        #   kinematic properties of contact point
        #   initial contact velocity
        if not self.initial_contact_velocity_calculated:
            self._dq0_n, self._dq0_t = self._contact_velocity(q)
            self.contact_model.set_dq0(self._dq0_n, self._dq0_t)
            self.initial_contact_velocity_calculated = True

        # if self.contact_detected:
        self._distance, self._delta = self._contact_geometry_GCS(q)
        # print "self._delta =", self._delta
        # print "self._dq0_n, self._dq0_t =", self._dq0_n, self._dq0_t

        #   current contact velocity at time t
        self._dq_n, self._dq_t = self._contact_velocity(q)
        # print "self._dq_n, self._dq_t =", self._dq_n, self._dq_t
        # time.sleep(100)
        #   current contact velocity at time t
        # self._dq_n, self._dq_t = self._contact_velocity(q)

        if self._type.lower() in self._types:#== "general" or self._type.lower() == "revolute clearance joint" or self._type.lower() == "contact sphere-sphere":#ECF-N
            self._solve_ECF_N(t, q, self._delta, self._dq_n, self._dq_t)#self._delta
        else:
            raise ValueError, "Contact type not correct!"


if __name__ == "__main__":
    a = RevoluteClearanceJoint(body_id_i=0, body_id_j=1, u_iP=np.array([0, 0]), u_jP=np.array([-1, 0]), type="rcj", properties_dict=[], parent=None)
    print "a =", a
    pprint(vars(a))
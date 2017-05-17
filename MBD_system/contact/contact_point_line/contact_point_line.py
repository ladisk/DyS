"""

created by: lskrinjar
date of creation: 13/03/2016
time of creation: 09:39
"""
import itertools
import time
import numpy as np
from OpenGL.GL import *


from MBD_system.contact.contact import Contact
from MBD_system.cad2cm_lcs import cad2cm_lcs
from MBD_system.contact.distance.distance_line_node import DistanceLineNode
from MBD_system.u_P_lcs2gcs import u_P_lcs2gcs
from MBD_system.t2n import t2n
from MBD_system.Ai_ui_P import Ai_ui_P_vector
from MBD_system.q2theta_i import q2theta_i
from MBD_system.q2R_i import q2R_i
from MBD_system.transform_cs import gcs2cm_lcs
from MBD_system.contact.evaluate_distance import evaluate_distance_2D
from MBD_system.n2t import n2t


class ContactPointLine(Contact):
    """
    classdocs
    """
    __id = itertools.count()

    def __init__(self, _type, body_id_i, body_id_j, u_iP, u_jP, u_jR, properties_dict=[], parent=None):
        """

        :param _type:
        :param body_id_i:   id of point body
        :param body_id_j:   id of line body
        :param u_iP:
        :param u_jP:
        :param u_jR:
        :param properties_dict:
        :param parent:
        :return:
        """
        #    number
        self.contact_id = self.__id.next()
        self._name = "Line_Point_Contact_"+str(self.contact_id)

        super(ContactPointLine, self).__init__(_type, body_id_i, body_id_j, name=self._name, properties_dict=properties_dict, parent=parent)


        #   type of contact
        self._contact_type = "contact point-line"

        #   predefine empty list of contact points
        self.u_P_LCS_list = []

        #   point on body i in body LCS
        self.u_iP_LCS = u_iP
        if self._parent is not None:
            if hasattr(self._parent._parent.bodies[body_id_i], "CM_CAD_LCS"):
                self.u_iP_LCS = cad2cm_lcs(u_iP, self._parent._parent.bodies[body_id_i].CM_CAD_LCS, 0)
        self.u_P_LCS_list.append(self.u_iP_LCS)

        #   line points on body j in body LCS
        self.u_jP_LCS = u_jP
        self.u_jR_LCS = u_jR
        if self._parent is not None:
            if hasattr(self._parent._parent.bodies[body_id_j], "CM_CAD_LCS"):
                self.u_jP_LCS = cad2cm_lcs(u_jP, self._parent._parent.bodies[body_id_j].CM_CAD_LCS, 0)
                self.u_jR_LCS = cad2cm_lcs(u_jR, self._parent._parent.bodies[body_id_j].CM_CAD_LCS, 0)

        self.u_P_LCS_list.append(self.u_jP_LCS)
        self.u_P_LCS_list.append(self.u_jR_LCS)

        #   append to line list
        self.u_P_line_LCS = []
        self.u_P_line_LCS.append(self.u_jP_LCS)
        self.u_P_line_LCS.append(self.u_jR_LCS)
        
        #    tangent
        self._t_LCS = self.u_P_line_LCS[1] - self.u_P_line_LCS[0]
        #   on line body j
        self._t_LCS_list[1] = self._t_LCS
        
        #    normal
        self._n_LCS = t2n(self._t_LCS)
        #   on line body j
        self._n_LCS_list[1] = self._n_LCS

        #   list of markers
        self.markers = self._create_markers()

    def _create_markers(self):
        """
        Function creates markers for point on body i and points that define line on body j
        :return:
        """
        markers = []
        #   create markers for points
        for i, u_P in enumerate(self.u_P_LCS_list):
            #   node coordinates
            node = np.array(np.append(u_P, self.z_dim), dtype="float32")

            if i == 0:
                body = self.body_list[0]
            if i == 1 or i == 2:
                body = self.body_list[1]

            #   create marker object
            marker = Marker(node, parent=body)

            #   append marker to list of body markers
            body.markers.append(marker)
            #   append marker to list of contact markers
            markers.append(marker)

        return markers

    def check_for_contact(self, q):
        """
        Function checks for contact for explicit point line contact object
        :param q:
        :return:
        """
        print "check_for_contact() ", self._type
        self._distance, self._delta = self._contact_geometry_GCS(q)
#         print "self._delta =", self._delta
        #   check sign
        self._sign_check = np.sign(self._delta * self._distance_solution_container[-1])

        if (self._sign_check == -1) and (self._delta <= 0):
            
            if abs(self._delta) <= self.distance_TOL:
                self._delta0 = self._delta
                print "self._delta0 =", self._delta0 
                self.contact_detected = True
                self.contact_distance_inside_tolerance = True
                self.status = 1
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

        return 0

    def contact_update(self, step, t, q):
        """
        Function updates status only for (general) contact
        :return:
            status - status of contact 0-no contact, 2-contact
        """
        self._step = step

        #   calculate relative contact velocity in normal direction at time t
        self._dq_n, self._dq_t = self._contact_velocity(q)

        #   calculate contact geometry in GCS to calculate contact distance
#         self._contact_geometry_GCS(q)

        #   calculate distance: node-edge
#         _distance, _inside = evaluate_distance_2D(self.node_GCS, self.edge_GCS[0], self._n_GCS, self._t_GCS)

        self._distance, self._delta = self._contact_geometry_GCS(q)
        # print "self._delta =", self._delta
        #   contact is present
        if self._delta < self._delta0:
            self.status = 1

        #   no contact
        else:
            self.status = 0
            self.contact_detected = False
            self.contact_distance_inside_tolerance = False
            # self.list_of_contact_force_objects_constructed = False
            self._contact_point_found = False
            self.initial_contact_velocity_calculated = False

            self.no_contact()

        # print "self.status - contact_update() =", self.status
        return self.status

    def _contact_geometry_GCS(self, q):
        """
        Function evaluates contact geometry in GCS and calculates distance from point to line
        :param q:
        :return:
        """
        #   point
        self.u_iP_GCS = u_P_lcs2gcs(self.u_iP_LCS, q, self.body_id_i)

        #   line
        self.u_P_line_GCS = []
        for u_P in self.u_P_line_LCS:
            u_P_GCS = u_P_lcs2gcs(u_P, q, self.body_id_j)
            self.u_P_line_GCS.append(u_P_GCS)

        #   normal
        self._n_GCS = Ai_ui_P_vector(self._n_LCS, q2theta_i(q, self.body_id_j))

        #   create distance object
        self._distance_obj = DistanceLineNode(self.u_iP_GCS, self.u_P_line_GCS[0], u_jR=self.u_P_line_GCS[1])

        #   distance and delta
        _distance = _delta = self._distance_obj._distance_sign

        return _distance, _delta

    def _get_contact_geometry_data(self, q):
        """
        
        contact point body i (point)
        contact point body j (line)
        """
        #    contact point on line
        #    GCS
        self.u_jP_GCS = self._distance_obj.contact_point_on_line()
        # print "self.u_jP_GCS =", self.u_jP_GCS
        # print "u_P_line_LCS =", self.u_P_line_LCS
        # time.sleep(100)
        
    def _contact_geometry_LCS(self, q):
        """
        
        """
        #    LCS
        self.u_jP_LCS = gcs2cm_lcs(self.u_jP_GCS, q2R_i(q, self.body_id_j), q2theta_i(q, self.body_id_j))
        
        #    list of contact points in LCS of each body
        self.u_P_LCS_list = [self.u_iP_LCS, self.u_jP_LCS]

        #   normal on body i - point
        self._n_LCS_list[0] = gcs2cm_lcs(-self._n_GCS, np.zeros(2), q2theta_i(q, self.body_id_i))

        #   tangent on body i - point
        self._t_LCS_list[0] = n2t(self._n_LCS_list[0])

    def _paint_GL_GCS(self, step=None):
        """
        Paint contact point in GCS
        :return:
        """
        if step is None:
            #   display contact point in GCS on each body in contact
            if self.contact_detected:
                for u_P_GCS, Fn, Ft in zip(self.u_P_GCS_list, self._Fn_list, self._Ft_list):
                    if Fn._visible or Ft._visible:
                        glVertex3f(u_P_GCS[0], u_P_GCS[1], self.z_dim)

            if self._distance_obj is not None:
                glColor3f(1, 0, 0)
                glVertex3f(self._distance_obj.r_iP[0], self._distance_obj.r_iP[1], self.z_dim)

                glColor3f(0, 0, 1)
                glVertex3f(self._distance_obj.r_jP[0], self._distance_obj.r_jP[1], self.z_dim)
                glVertex3f(self._distance_obj.r_jR[0], self._distance_obj.r_jR[1], self.z_dim)

            # for distance_obj in self._distance_list:
            #     distance_obj._paint_GL()
        else:
            None
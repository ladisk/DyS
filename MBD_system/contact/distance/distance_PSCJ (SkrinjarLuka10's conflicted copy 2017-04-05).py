"""

created by: lskrinjar
date of creation: 17/01/2017
time of creation: 20:33
"""
import time
import itertools
import numpy as np
from matplotlib import pyplot as plt
from pprint import pprint


from MBD_system.A import A_matrix
from MBD_system.transform_cs import gcs2cm_lcs, cm_lcs2gcs
from MBD_system.q2R_i import q2R_i
from MBD_system.q2theta_i import q2theta_i
from MBD_system.q2dR_i import q2dR_i
from MBD_system.q2dtheta_i import q2dtheta_i
from MBD_system.n2t import n2t
from MBD_system.dr_contact_point_uP import dr_contact_point_uP
from MBD_system.contact.distance.distance import Distance
from MBD_system.contact.distance.distance_line_node import DistanceLineNode


class DistancePSCJ(object):
    """
    classdocs
    Class object of distance line to node for PSCJ (Pin-Slot Clearance Joint)
    """
    _id = itertools.count(0)

    def __init__(self, r_iP, r_iR, r_jP, normal=None, tangent=None, R0_j=0., h0_iP=0., c=0., L=0., edge_num=None, parent=None):
        """
        :param r_iP:        line point in GCS on body i - slot
        :param r_iR:        line point in GCS on body i - slot
        :param r_jP:        free a point in GCS on body j - pin
        :param normal:      normal in GCS
        :param tangent:     tangent in GCS
        """
        #   id
        self.id = self._id.next()

        #   parent
        self._parent = parent

        #    color
        if hasattr(self._parent, "color"):
            self.color = self._parent.color
        else:
            self.color = np.random.rand(3)

        #   body ids
        #   i - free point of pin body id
        #   j - line body of slot body id
        if self._parent is not None:
            self.body_id_i = self._parent.body_id_i
            self.body_id_j = self._parent.body_id_j
        else:
            self.body_id_i = None
            self.body_id_j = None

        self.body_id_list = [self.body_id_i, self.body_id_j]

        #   coordinates in GCS
        self.r_iP = r_iP
        self.r_iR = r_iR
        self.r_jP = r_jP

        #   index of which contact data from list is used at contact point
        self.index = None
        # print "self.r_iP =", self.r_iP
        # print "self.r_iR =", self.r_iR
        # print "self.r_jP =", self.r_jP

        self.r_iPR_list = [self.r_iP, self.r_iR]
        self.r_jPiL_list = [np.zeros(2, dtype=float), np.zeros(2, dtype=float)]

        self.sign_list = [+1, -1]

        self.sign_cylindrical_section = [-1, +1]

        #   radius of a pin body j
        if self._parent is not None:
            self.R0_j = self._parent.R0_j
            self.c = self._parent._radial_clearance
            self.L = L
            self.h0_iP = self._parent.h0_iP

        else:
            self.R0_j = R0_j
            self.c = c
            self.L = L
            self.h0_iP = h0_iP

        self.R0_i = self.h0_iP / 2.

        #   coordinates in LCS
        # self.u_iP = u_iP
        # self.u_jP = u_jP
        # self.u_jR = u_jR

        #   edge number
        self.edge_num = edge_num

        #   normal, tangent in LCS
        # self.tangent_LCS = tangent_LCS
        # self.normal_LCS = normal_LCS
        self._contact_point_found = False

        #   edge vector
        self.edge = r_iR - r_iP

        #   tangent
        if tangent is None and normal is None and self.edge is not None:
            self.tangent = self.edge / np.linalg.norm(self.edge, ord=2)
        elif tangent is None and normal is not None:
            self.tangent = np.dot(A_matrix(np.pi / 2), normal)
        else:
            self.tangent = tangent

        #   normal
        if normal is None and self.tangent is not None:
            self.normal = np.array([self.tangent[1], -self.tangent[0]])
        else:
            self.normal = normal

        #   z dimension
        self.z_dim = 0.

        #   list of distance vector values
        self._distance_vector = [np.zeros(2, dtype=float),
                                 np.zeros(2, dtype=float),
                                 np.zeros(2, dtype=float),
                                 np.zeros(2, dtype=float)]

        #   list of position options of contact point
        self.contact_point_position_list = ["AiBi",
                                            "CiDi",
                                            "iPjP",
                                            "iRjP"]

        self.sections = ["jPjR",
                         "jP",
                         "jR"]

        self.pin_in_section = None

        self.phi_list = np.zeros(2, dtype=float)

        self.n_list = 4 * [None]
        self.t_list = 4 * [None]

        #   distance vector
        self._distance_vector = self._evaluate_distance_vector()

        #   evaluate pin position in a slot
        self.pin_in_section = self.evaluate_pin_in_section()

        #   distance and distance sign
        self._evaluate_distance()

    # def _evaluate_angles(self):
    #     """
    #     Function evaluates angles fi1, fi2 between:
    #     fi1 - edge vector and distance vector
    #     fi2 - edge vector and vector rPi- rRj
    #     :return:
    #     """
    #     fi1 = np.arccos(np.dot(self.tangent, self._distance_vector) / np.linalg.norm(self._distance_vector, ord=2))
    #     self._distance_vector_2 = self._evaluate_distance_vector(self.r_jP, self.r_iR)
    #     fi2 = np.arccos(np.dot(-self.tangent, self._distance_vector_2) / np.linalg.norm(self._distance_vector_2, ord=2))
    #
    #     return fi1, fi2

    def _evaluate_distance_vector(self):
        """
        Evaluate distance vector from
        :return:
        """
        # print "_evaluate_distance_vector()@",self
        distance_vector = [np.zeros(2, dtype=float),
                           np.zeros(2, dtype=float),
                           np.zeros(2, dtype=float),
                           np.zeros(2, dtype=float)]

        #   distance vector at flat section
        for i, (r_iL, sign) in enumerate(zip(self.r_iPR_list, self.sign_list)):
            t = sign * self.tangent
            # print "t =", t
            self.t_list[i] = t
            # print "i =", i
            # print "t =", self.t_list[i]
            self.n_list[i] = np.array([-t[1], t[0]])
            # print "n =", self.n_list[i]
            # print "r_iL =", r_iL
            # print "self.r_iP - r_iL =",

            r_jPiL = self.r_jP - r_iL
            self.r_jPiL_list[i] = r_jPiL
            # print "r_jPiL =", r_jPiL
            distance_vector[i] = np.dot(r_jPiL, self.n_list[i]) * self.n_list[i]
            # print "d =", np.linalg.norm(distance_vector[i])
            # print "TEST =", np.dot(distance_vector[i], n)

            # self.phi_list[i] = np.arctan2(r_jPiL[1] - self.t_list[i][1], r_jPiL[0] - self.t_list[i][0])
            # self.phi_list[i] = np.arccos(np.dot(r_jPiL, self.tangent)/np.linalg.norm(r_jPiL))
            # self.phi_list[i] = np.arctan2(r_jPiL[1], r_jPiL[0]) - np.arctan2(t[1], t[0])

        # print "distance_vector at flat section ="
        # print distance_vector

        #   distance vector at cylindrical section
        for i, r_iL in enumerate(self.r_iPR_list):
            d_vec = self.r_jP - r_iL
            distance_vector[i + 2] = d_vec
        # self.t_list[i] = t
        # # print "i =", i
        # # print "t =", self.t_list[i]
            n_i = d_vec / np.linalg.norm(d_vec)
            self.n_list[i + 2] = n_i

            t_i = np.array([-n_i[1], n_i[0]])
            self.t_list[i + 2] = t_i
        #   evaluate phi at each cylindrical section
        self.phi_list = self._evaluate_phi(distance_vector)

        return distance_vector

    def evaluate_distance_vector(self):
        """

        :return:
        """

    def evaluate_pin_in_section(self):
        """

        :return:
        """
        pin_in_section = ""

        if (np.pi / 2. < self.phi_list[0] < 3. * np.pi / 2.) and (np.pi / 2. < self.phi_list[1] < 3. * np.pi / 2.):
            pin_in_section = "iPiR"

        else:
            for phi, section in zip(self.phi_list, ["iP", "iR"]):
                if -np.pi / 2. <= phi <= np.pi / 2.:
                    pin_in_section = section

        return pin_in_section

    def _evaluate_distance(self):
        """
        Variable distance (=delta) is negative when in contact!
        :return:
        """
        #   predefined empty list
        self._distance = 4 * [None]
        self._contact = 4 * [False]

        #   deformation at flat section
        for i in range(0, 2):
            n_i = self.n_list[i]
            # print "n_i =", n_i
            d_vec = self._distance_vector[i] - n_i * self.R0_i
            d = np.linalg.norm(d_vec)
            # print "d =", d
            self._distance[i] = d - self.R0_j
            # print "self._distance[i] =", self._distance[i]
            # print "self.n_list[i] =", self.n_list[i]
            # print "1 =", np.dot(self._distance_vector[i], self.n_list[i]) > 0.
            # print "2 =", self._distance[i] > 0., self._distance[i]
            if (np.pi / 2. < self.phi_list[0] < 3. * np.pi / 2.) and (np.pi / 2. < self.phi_list[1] < 3. * np.pi / 2.):
                self.pin_in_section = "iPiR"
                if (np.dot(self._distance_vector[i], self.n_list[i]) > 0.) and (self._distance[i] > 0.):
                    self._contact[i] = True

        #   deformation at cylindrical section
        for i, phi in enumerate(self.phi_list):
            d = np.linalg.norm(self._distance_vector[i + 2])
            self._distance[i + 2] = self.c - d

            if -np.pi / 2. <= phi <= np.pi / 2.:
                if i == 0:
                    self.pin_in_section = "iP"
                else:
                    self.pin_in_section = "iR"

                if d > 0.:
                    self._contact[i + 2] = True

        # print "self._distance ="
        # print self._distance
        #
        # print "self._contact ="
        # print self._contact


        # for i, d in enumerate(self._distance):
        #     print "i =", i, "d =", d
        # #   distance projection
        # _distance_vector_projection = self.normal * np.dot(self._distance_vector, self.normal)
        #
        # #   distance
        # self._distance = np.linalg.norm(_distance_vector_projection, ord=2)
        # for i, d in enumerate(self._distance):
        #     if d < 0.:
        #         self._inside[i] = True

        #    distance sign
        # self._inside = False
        # self._distance_sign = self._distance

        # if (np.sign(np.cross(self.tangent, self._distance_vector)) >= 0) and (self.r_iR is not None):
            # self.fi1, self.fi2 = self._evaluate_angles()
            # print 'self.fi1=', np.rad2deg(self.fi1), 'self.fi2=', np.rad2deg(self.fi2), self.fi1 <= np.pi/2. and self.fi2 <= np.pi/2.
            # if self.fi1 <= np.pi/2. and self.fi2 <= np.pi/2.:
        # self._inside = True
        # self._distance_sign = -self._distance
        # print '_evaluate_distance=', self._distance_sign

    def _evaluate_phi(self, distance_vector):
        """

        :return:
        """
        # print "_evaluate_phi()"
        # print "distance_vector =", distance_vector
        phi_list = []
        for i, sign in enumerate(self.sign_cylindrical_section):
            # print "----------------------------------"
            r_jPiL = distance_vector[i + 2]
            # print "r_jPiL =", r_jPiL
            t_i = sign * self.tangent
            # print "t_i =", t_i
            # print "r_jPiL =", r_jPiL
            # print "fi_t =", np.rad2deg(np.arctan2(t[1], t[0]))
            # print "fi_r =", np.rad2deg(np.arctan2(r_jPiL[1], r_jPiL[0]))
            # print "test =", np.rad2deg(np.arccos(np.dot(r_jPiL, t)/np.linalg.norm(r_jPiL)))
            # print "test =", np.rad2deg(np.arctan2(r_jPiL[1], r_jPiL[0]) - np.arctan2(t[1], t[0]))
            # angles[i] = np.arctan2(t[1] - r_jPiL[1], t[0] - r_jPiL[0])
            # phi_i = np.arctan2(r_jPiL[1], r_jPiL[0]) - np.arctan2(t_i[1], t_i[0])
            phi_i = np.arccos(np.dot(r_jPiL, t_i)/np.linalg.norm(r_jPiL))
            # print "phi_i =", np.rad2deg(phi_i)
            phi_list.append(phi_i)

        return phi_list

    def _evaluate_CP_GCS(self):
        """

        :return:
        """
        self.CP = self.r_iP + np.dot(self.r_iPjP, self.tangent) * self.tangent

    def contact_point_on_line(self):
        """
        Function calculates contact point - CP
        :return:
        """
        self._evaluate_CP_GCS()
        return self.CP

    def geometry_LCS(self):
        """
        Method transforms all 3 vectors (free node + edge nodes) from GCS to LCS
        :param q:
        :param point_body_id:
        :param edge_body_id:
        :return:
        """
        for r, u, body_id in zip([self.r_iP, self.r_jP, self.r_jR], ["u_iP", "u_jP", "u_jR"], [self.body_id_i, self.body_id_j, self.body_id_j]):
            _u = self._parent._contact_geometry_LCS(body_id, r)

            setattr(self, u, _u)

        #   normal
        self.normal_LCS = self._parent._contact_normal_LCS(self.body_id_j, self.normal)

        #   tangent
        self.tangent_LCS = self._parent._contact_tangent_LCS(self.body_id_j, self.tangent)

    def update_contact_geometry_GCS(self, q, u_iP=None, u_jP=None, normal_LCS=None, u_jR=None, tangent_LCS=None):
        """
        Method updates - calculates contact geometry coordinates of all 3 vector (free node + edge nodes) from LCS to GCS
        :param q_i: coordinates of a point body
        :param q_j: coordinates of a edge body
        :return:
        """
        if self.u_iP is None and u_iP is not None:
            self.u_iP = u_iP

        if self.u_jP is None and u_jP is not None:
            self.u_jP = u_jP

        if self.u_jR is None and u_jR is not None:
            self.u_jR = u_jR

        if self.normal_LCS is None and normal_LCS is not None:
            self.normal_LCS = normal_LCS

        if self.tangent_LCS is None and tangent_LCS is not None:
            self.tangent_LCS = tangent_LCS

        #   to GCS
        for u, r, body_id in zip([self.u_iP, self.u_jP, self.u_jR], ["r_iP", "r_jP", "r_jR"], [self.body_id_i, self.body_id_j, self.body_id_j]):
            R = q2R_i(q, body_id)
            theta = q2theta_i(q, body_id)

            #   point coordinate in LCS
            _r = cm_lcs2gcs(u, R, theta)

            #   set object attribute
            setattr(self, r, _r)

        theta = q2theta_i(q, self.body_id_j)

        #   normal in LCS
        self.normal = cm_lcs2gcs(self.normal_LCS, np.zeros(2), theta)

        #   tangent in LCS
        self.tangent = n2t(self.normal)

        #   distance
        self._distance_vector = self._evaluate_distance_vector(self.r_iP, self.r_jP)

        #   distance and distance sign
        self._evaluate_distance()

    def inside_penetration_depth(self, max_penetration_depth):
        """
        Function checks if distance between free node and edge is inside maximum penetration depth
        :param max_penetration_depth:
        :return:
        """
        bc = self.r_iP - self.r_jP
        ba = self.r_iP - self.r_jP

        #   dot product
        dotbac = np.dot(bc, ba)

        #   theta angle
        self.theta = np.rad2deg(np.arccos(dotbac / (np.linalg.norm(bc, ord=2) * np.linalg.norm(ba, ord=2))))

        if (self.theta <= 180.) and (self._distance_sign < 0.) and (abs(self._distance_sign) <= max_penetration_depth):
            return True

        return False

    def contact_geometry_GCS(self):
        """

        :return:
        """
        # print "contact_geometry_GCS()@",__name__
        # print "self.pin_in_section =", self.pin_in_section
        # print "self._distance[0:2] =", self._distance[0:2]
        if self.pin_in_section == "iPiR":
            self.index = self._distance[0:2].index(min(self._distance[0:2]))
            # print "self.index =", self.index
            distance = np.linalg.norm(self._distance_vector[self.index])
            delta = self._distance[self.index]

            t_GCS = self.tangent * self.sign_list[self.index]
            n_GCS = np.array([t_GCS[1], -t_GCS[0]])

        if self.pin_in_section in ["iP", "iR"]:
            if self.pin_in_section == "iP":
                self.index = 2
            if self.pin_in_section == "iR":
                self.index = 3

            distance = self._distance_vector[self.index]
            delta = self._distance[self.index]

            n_GCS = distance / np.linalg.norm(distance)
            t_GCS = np.array([-n_GCS[1], n_GCS[0]])

        return distance, delta, n_GCS, t_GCS


if __name__ == "__main__":
    r_iP = np.array([0., 0.])
    r_iR = np.array([1E-3, 0.])
    r_jP = np.array([0.5E-3, -.2E-3])

    # r_iP = np.array([-0.000833658094, -0.000216038301], dtype=float)
    # r_iR = np.array([8.482413250022e-05, 2.663517625479e-05], dtype=float)
    # r_jP = np.array([9.999999747379e-06, -4.199999966659e-05], dtype=float)

    L = 0.95E-3
    R0j = 0.75E-3#0.1
    hi = 1.6E-3#0.3
    R0i = hi / 2.
    c = R0i - R0j
    d_PSCJ = DistancePSCJ(r_iP, r_iR, r_jP, R0_j=R0j, h0_iP=hi, c=c, L=L)

    fig = plt.figure(num=1, figsize=(6, 5), dpi=72, facecolor='w', edgecolor='g')
    ax = plt.subplot(111, aspect="equal")

    #   pin - body j
    plt.plot(r_jP[0], r_jP[1], color='r')
    circle_j = plt.Circle((r_jP[0], r_jP[1]), R0j, color='r', fill=False)
    ax.add_artist(circle_j)

    #   slot - body i
    circle_iP = plt.Circle((r_iP[0], r_iP[1]), R0i, color='b', fill=False)
    ax.add_artist(circle_iP)
    circle_iR = plt.Circle((r_iR[0], r_iR[1]), R0i, color='b', fill=False)
    ax.add_artist(circle_iR)
    #   tangent - center line
    plt.plot([r_iP[0], r_iR[0]], [r_iP[1], r_iR[1]], color="black")

    #   plot points
    for r in [r_iP, r_iR, r_jP]:
        plt.plot(r[0], r[1], marker="o", color="black")

    plt.xlim([min(r_iP[0], r_iR[0]) - 2 * R0i, max(r_iP[0], r_iR[0]) + 2 * R0i])
    plt.ylim([min(r_iP[1], r_iR[1]) - 2 * R0i, max(r_iP[1], r_iR[1]) + 2 * R0i])

    plt.show()
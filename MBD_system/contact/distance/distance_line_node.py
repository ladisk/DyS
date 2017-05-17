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


class DistanceLineNode(object):#Distance
    """
    classdocs
    """
    _id = itertools.count(0)

    def __init__(self, r_iP, r_jP, normal=None, r_jR=None, tangent=None, u_iP=None, u_jP=None, normal_LCS=None, u_jR=None, tangent_LCS=None, body_id_i=None, body_id_j=None, edge_num=None, theta_jP=None, theta_jR=None, parent=None):
        """
        :param r_iP:        free a point in GCS on body i
        :param r_jP:        line point in GCS on body j
        :param normal:      normal in GCS
        :param r_jR:        line point in GCS on body j
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
        #   i - free point of body id
        #   j - line body of body id
        self.body_id_i = body_id_i
        self.body_id_j = body_id_j
        self.body_id_list = [self.body_id_i, self.body_id_j]

        #   coordinates in GCS
        self.r_iP = r_iP
        self.r_jP = r_jP
        self.r_jR = r_jR

        #   coordinates in LCS
        self.u_iP = u_iP
        self.u_jP = u_jP
        self.u_jR = u_jR

        #   theta angles
        self.theta_jP = theta_jP
        self.theta_jR = theta_jR
        if self.theta_jP is None:
            self.theta_jP = 180.

        if self.theta_jR is None:
             self.theta_jR = 180.

        #   max penetration depth
        self.max_penetration_depth = 5E-4
        if self._parent is not None:
            if hasattr(self._parent._parent._parent._parent, "bodies"):
                self.max_penetration_depth = self._parent._parent._parent._parent.bodies[self.body_id_j].max_penetration_depth

        #   edge number
        self.edge_num = edge_num

        #   normal, tangent in LCS
        self.tangent_LCS = tangent_LCS
        self.normal_LCS = normal_LCS
        self._contact_point_found = False

        #   contact status
        self._inside = False

        #   edge vector
        if r_jR is not None:
            self.edge = r_jR - r_jP
        else:
            self.edge = None

        #   tangent
        if tangent is None and normal is None and self.edge is not None:
            self.tangent = self.edge / np.linalg.norm(self.edge, ord=2)
        elif tangent is None and normal is not None:
            self.tangent = np.dot(A_matrix(np.pi/2), normal)
        else:
            self.tangent = tangent

        #   normal
        if normal is None and self.tangent is not None:
            self.normal = np.array([self.tangent[1], -self.tangent[0]])
        else:
            self.normal = normal

        #   z dimension
        self.z_dim = 0.

        if all(self.r_iP == np.zeros(2)) and all(self.r_jP == np.zeros(2)):
            pass
        else:
            #   distance vector
            self._distance_vector = self._evaluate_distance_vector(self.r_iP, self.r_jP)

            #   distance and distance sign
            self._evaluate_distance()

    def _evaluate_angles(self):
        """
        Function evaluates angles fi1, fi2 between:
        fi1 - edge vector and distance vector
        fi2 - edge vector and vector rPi- rRj
        :return:
        """
        print "t =", self.tangent, id(self)
        fi1 = np.arccos(np.dot(self.tangent, self._distance_vector) / np.linalg.norm(self._distance_vector, ord=2))
        self._distance_vector_2 = self._evaluate_distance_vector(self.r_iP, self.r_jR)
        fi2 = np.arccos(np.dot(-self.tangent, self._distance_vector_2) / np.linalg.norm(self._distance_vector_2, ord=2))

        return fi1, fi2

    def _evaluate_distance_vector(self, r_iP, r_jP):
        """

        :return:
        """
        _distance_vector = r_iP - r_jP

        return _distance_vector

    def _evaluate_distance(self):
        """

        :return:
        """

        r_BA = self.r_jR - self.r_jP
        r_CA = self.r_iP - self.r_jP

        #   distance projection
        _distance_vector_projection = self.normal * np.dot(self._distance_vector, self.normal)

        #   distance
        self._distance = np.linalg.norm(_distance_vector_projection, ord=2)

        #    distance sign
        self._distance_sign = self._distance

        # print 'r_iP=', self.r_iP
        # print 'r_jP=', self.r_jP
        # print 'r_jR=', self.r_jR
        # print '.................'

        #   State no.:1
        # r_BA = self.r_jR - self.r_jP
        # r_CA = self.r_iP - self.r_jP

        BACA = np.cross(r_BA, r_CA)
        # BACA = np.cross(r_CA, r_BA)
        # print 'State_1=', BACA >=0., BACA, self.normal

        #   State no.:2
        r_AB = self.r_jP - self.r_jR
        r_AC = self.r_jP - self.r_iP

        dotABAC = np.dot(r_AB, r_AC)
        thetaA = np.rad2deg(np.arccos(dotABAC / (np.linalg.norm(r_AB, ord=2) * np.linalg.norm(r_AC, ord=2))))
        # print 'thetaA=', thetaA

        #   State no.:3
        r_AB = self.r_jP - self.r_jR
        r_CB = self.r_iP - self.r_jR
        dotABCB = np.dot(r_AB, r_CB)
        thetaB = np.rad2deg(np.arccos(dotABCB / (np.linalg.norm(r_AB, ord=2) * np.linalg.norm(r_CB, ord=2))))
        # cond3 = self.thetaB <= (self.theta_jR / 2.)
        # print 'r_BA=', r_BA, 'r_CA=', r_CA
        # print 'cond3=', self.cond3, self.thetaB, self.theta_jR / 2.

        #   State no.:4

        #   penetration state
        cond1 = BACA >= 0.

        # print "self._parent =", self._parent
        cond2 = thetaA <= (self.theta_jP / 2.)
        # print "cond2 =", cond2, thetaA, self.theta_jP / 2.
        cond3 = thetaB <= (self.theta_jR / 2.)
        # print "cond3 =", cond3
        cond4 = abs(self._distance_sign) <= self.max_penetration_depth
        # print "cond4 =", cond4
        if cond1 and cond2 and cond3 and cond4:
            self._inside = True
            self._distance_sign = -self._distance
        else:
            self._inside = False

    def inside(self):
        """

        :return:
        """
        self._evaluate_distance()
        # print "self._parent =", self._parent
        return self._inside

    def _evaluate_CP_GCS(self):
        """

        :return:
        """
        self.CP = self.r_jP + np.dot(self._distance_vector, self.tangent) * self.tangent #self.r_iP +

    def contact_point_on_line(self):
        """
        Function calculates contact point - CP
        :return:
        """
        self._evaluate_CP_GCS()
        return self.CP

    def geometry_LCS(self, q=None):
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

        if q is not None:
            #   update global coordinates
            self.update_contact_geometry_GCS(q)

    def update_contact_geometry_GCS(self, q, u_iP=None, u_jP=None, normal_LCS=None, u_jR=None, tangent_LCS=None):
        """
        Method updates - calculates contact geometry coordinates of all 3 vector (free node + edge nodes) from LCS to GCS
        :param q_i: coordinates of a point body
        :param q_j: coordinates of a edge body
        :return:
        """
        # print 'q=', q
        if self.u_iP is None and u_iP is not None:
            self.u_iP = u_iP
        # print 'self.u_iP=', self.u_iP

        if self.u_jP is None and u_jP is not None:
            self.u_jP = u_jP
        # print 'self.u_jP=', self.u_jP

        if self.u_jR is None and u_jR is not None:
            self.u_jR = u_jR
        # print 'self.u_jR=', self.u_jR

        if self.normal_LCS is None and normal_LCS is not None:
            self.normal_LCS = normal_LCS

        if self.tangent_LCS is None and tangent_LCS is not None:
            self.tangent_LCS = tangent_LCS

        #   to GCS
        for u, r, body_id in zip([self.u_iP, self.u_jP, self.u_jR], ["r_iP", "r_jP", "r_jR"], [self.body_id_i, self.body_id_j, self.body_id_j]):
            R = q2R_i(q, body_id)
            theta = q2theta_i(q, body_id)

            #   point coordinate in LCS
            # print 'u=',u
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

    def plot(self, show=False):
        """

        Returns:
        """
        fig = plt.figure(num=1, figsize=(6, 5), dpi=200, facecolor='w', edgecolor='k')
        ax = plt.subplot(111, aspect="equal")
        ax.ticklabel_format(style='sci', axis='both')

        #   free node
        plt.plot(self.r_iP[0], self.r_iP[1], 'D')
        plt.text(self.r_iP[0], self.r_iP[1], str(self.id))
        plt.text(self.r_iP[0]*1.1, self.r_iP[1]*1.1, "[" + str(round(self.r_iP[0]*1E+3, 3)) + "," + str(round(self.r_iP[0]*1E+3, 3)) + "]", fontsize=8)

        #   egde
        plt.plot([self.r_jP[0], self.r_jR[0]], [self.r_jP[1], self.r_jR[1]])

        #   AABB
        #   TODO - check
        # print "self.overlap_pair =", self.overlap_pair
        if hasattr(self, "overlap_pair"):
            if hasattr(self.overlap_pair, "_AABB_list"):
                AABB = self.overlap_pair._AABB_list[1]
                AABB.plot_2D(_ax=ax)

        if show:
            fig.show()

    def save_plot(self, step=None):
        """

        Returns:
        """
        fig = plt.figure(num=1, figsize=(6, 4), dpi=100, facecolor='w', edgecolor='k')
        ax = plt.subplot(111, aspect="equal")
        if step is None:
            filename = str(self.id)+".png"
        else:
            filename = "step_"+str(step)+".png"

        self.plot()
        plt.xlim([-0.025, +0.009])
        plt.ylim([-0.009, +0.009])
        plt.savefig(filename)
        print "Plot saved to file: ", filename
        plt.clf()


if __name__ == "__main__":
#     #   free node
#     r_iP = np.array([1.5, -1.1346])
#     #   line node
#     r_jP = np.array([0, 0])
#     normal = np.array([0, 1])
#     d = DistanceLineNode(r_iP, r_jP, normal)
#     print d._distance

    fig = plt.figure(num=1, figsize=(6, 4), dpi=100, facecolor='w', edgecolor='k')
    ax = plt.subplot(111, aspect="equal")
#     #   normal
#     plt.quiver([normal[0], normal[0]], [0, normal[1]], linestyle="-", color=(0, 0, 0))
#     #   line node
#     plt.plot(r_jP[0], r_jP[1], marker="x", color=(0, 0, 0))
#     #   free node
#     plt.plot(r_iP[0], r_iP[1], marker="x", color=(1, 0, 0), markersize=5)
#     #   distance
#     plt.plot([r_iP[0], r_iP[0]], [0, d._distance], linestyle="-", color=(0, 1, 0), markersize=5)
# 
#     plt.show()

    #   free node
    r_iP = np.array([4.316930274e-08,  5.235721891e-05]) #OK
    # r_iP = np.array([0.011439999485,  0.006928203106]) #NOK
    #   line node
    r_jP = np.array([0.000000000e+00,  3.000000000e-03])
    r_jR = np.array([0.000000000e+00, -3.000000000e-03])

    #   free node
    # r_iP = np.array([-0.007592, -0.011551])
    # # r_iP = np.array([.4, .6])
    # #   line node
    # r_jR = np.array([-0.007, -0.01])
    # r_jP = np.array([-0.009, -0.015])
    
    d = DistanceLineNode(r_iP, r_jP, r_jR=r_jR)
    d.normal = np.array([-1., 0.])
    print d._distance
    
    plt.plot([r_jP[0], r_jR[0]], [r_jP[1], r_jR[1]])
    plt.plot(r_iP[0], r_iP[1], "o")
    
    proj = d.contact_point_on_line()
    print "point on line =", proj
    pprint(vars(d))
    print d._evaluate_angles()
    # plt.plot(proj[0], proj[1], "s")
    # plt.show()




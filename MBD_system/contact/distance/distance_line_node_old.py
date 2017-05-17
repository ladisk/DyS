__author__ = 'lskrinjar'
import time
import itertools
import numpy as np
from matplotlib import pyplot as plt
from OpenGL.GL import *
from OpenGL.GLU import *
from pprint import pprint


from MBD_system.A import A_matrix
from MBD_system.transform_cs import gcs2cm_lcs, cm_lcs2gcs
from MBD_system.q2R_i import q2R_i
from MBD_system.q2theta_i import q2theta_i
from MBD_system.q2dR_i import q2dR_i
from MBD_system.q2dtheta_i import q2dtheta_i
from MBD_system.n2t import n2t
from MBD_system.t2n import t2n
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
            self.normal = t2n(self.tangent)
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
        # print "t =", self.tangent, id(self)
        fi1 = np.arccos(np.dot(self.tangent, self._distance_vector) / np.linalg.norm(self._distance_vector, ord=2))
        self._distance_vector_2 = self._evaluate_distance_vector(self.r_iP, self.r_jR)
        fi2 = np.arccos(np.dot(-self.tangent, self._distance_vector_2) / np.linalg.norm(self._distance_vector_2, ord=2))

        return fi1, fi2

    def _evaluate_distance_vector(self, r_iP, r_jP):
        """

        :return:
        """
        # Problem, ker je pred kontaktom r_jp robna tocka pol ob zacetku kontakta pa postane kontaktna tocka
        _distance_vector = r_iP - r_jP
        # print 'r_iP=', r_iP, 'r_jP=', r_jP
        # print 'DV=',_distance_vector

        return _distance_vector

    def _evaluate_distance(self):
        """

        :return:
        """
        #   distance projection
        _distance_vector_projection = self.normal * np.dot(self._distance_vector, self.normal)
        # print 'normal=', self.normal
        # print 'distance_vector=', self._distance_vector
        # print 'DVP=',_distance_vector_projection

        #   distance
        self._distance = np.linalg.norm(_distance_vector_projection, ord=2)
        # print 'DIST=', self._distance

        #    distance sign
        self._distance_sign = self._distance

        #   State no.:1
        r_BA = self.r_jR - self.r_jP
        r_CA = self.r_iP - self.r_jP

        self.crossBACA = np.cross(r_BA, r_CA)
        self.cond1 = self.crossBACA >= 0.

        #   State no.:2
        r_AB = self.r_jP - self.r_jR
        r_AC = self.r_jP - self.r_iP

        dotABAC = np.dot(r_AB, r_AC)
        self.thetaA = np.rad2deg(np.arccos(dotABAC / (np.linalg.norm(r_AB, ord=2) * np.linalg.norm(r_AC, ord=2))))
        self.cond2 = self.thetaA <= (self.theta_jP / 2.)
        # print 'cond2=', self.cond2, self.thetaA, self.theta_jP / 2.

        #   State no.:3
        r_AB = self.r_jP - self.r_jR
        r_CB = self.r_iP - self.r_jR
        dotABCB = np.dot(r_AB, r_CB)
        self.thetaB = np.rad2deg(np.arccos(dotABCB / (np.linalg.norm(r_AB, ord=2) * np.linalg.norm(r_CB, ord=2))))
        self.cond3 = self.thetaB <= (self.theta_jR / 2.)
        # print 'r_BA=', r_BA, 'r_CA=', r_CA
        # print 'cond3=', self.cond3, self.thetaB, self.theta_jR / 2.

        #   State no.:4
        self.cond4 = abs(self._distance_sign) <= self.max_penetration_depth
        # print 'cond4=', self.cond4, self._distance_sign, self.max_penetration_depth

        #   Penetration condition
        if self.cond1 and self.cond2 and self.cond3 and self.cond4:
            self._inside = True
            # print 'inside = TRUE', self.r_iP, self.r_jP, self.r_jR
            self._distance_sign = -self._distance
        else:
            self._inside = False

        # if self._inside == True:
            # print 'INS=', self._inside
            # print 'DVP=', _distance_vector_projection
            # print 'delta=', self._distance_sign

    def inside(self):
        """

        :return:
        """
        self._evaluate_distance()

        return self._inside

    def _evaluate_CP_GCS(self):
        """

        :return:
        """
        self.CP = self.r_jP + np.dot(self._distance_vector, self.tangent) * self.tangent
        # print 'r_jP=', self.r_jP
        # print 'self._distance_vector=',self._distance_vector
        # print 'self.tangent=', self.tangent
        # print 'CP @ DLN=', self.CP
        # print '...................'

    def contact_point_on_line(self):
        """
        Function calculates contact point - CP
        :return:
        """
        self._evaluate_CP_GCS()

        return self.CP

    def geometry_LCS(self, q):
        """
        Method transforms all 3 vectors (free node + edge nodes) from GCS to LCS
        :param q:
        :param point_body_id:
        :param edge_body_id:
        :return:
        """
        # Checked, working OK!
        for r, u, body_id in zip([self.r_iP, self.r_jP, self.r_jR], ["u_iP", "u_jP", "u_jR"], [self.body_id_i, self.body_id_j, self.body_id_j]):
            _u = self._parent._contact_geometry_LCS(body_id, r)

            setattr(self, u, _u)

        #   normal
        self.normal_LCS = self._parent._contact_normal_LCS(self.body_id_j, self.normal)
        # print "self.normal_LCS =", self.normal_LCS

        #   tangent
        self.tangent_LCS = self._parent._contact_tangent_LCS(self.body_id_j, self.tangent)
        # print "self.tangent_LCS =", self.tangent_LCS

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

            #   point coordinate in GCS
            _r = cm_lcs2gcs(u, R, theta)
            # print "u =", u, "r =", _r
            #   set object attribute
            setattr(self, r, _r)

        theta = q2theta_i(q, self.body_id_j)

        #   normal in GCS
        self.normal = cm_lcs2gcs(self.normal_LCS, np.zeros(2), theta)
        # print "self.normal =", self.normal
        #   tangent in GCS
        self.tangent = n2t(self.normal)
        # print "self.tangent =", self.tangent
        #   distance
        self._distance_vector = self._evaluate_distance_vector(self.r_iP, self.r_jP)
        # print "self._distance_vector =", self._distance_vector
        #   distance and distance sign
        self._evaluate_distance()
        # print "distance_sign =", self._distance_sign


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
    # test 1
    #   free node - brush ; theta=0,0,0
    r_iP = np.array([4.565690181e-09,  0.000000000e+00]) #OK
    r_ix= np.linspace(0,1.e-8,1000)

    #   line node - slip-ring; theta=0,0,1
    r_jP = np.array([-5.235721891e-05,  2.999543085e-03])
    r_jR = np.array([5.235721891e-05, -2.999543085e-03])

    # test 20
    # free node - brush ; theta=0,0,0
    r_iP = np.array([4.565690181e-09,  0.000000000e+00])  # OK
    r_ix = np.linspace(0, 1.e-6, 1000)

    #   line node - slip-ring; theta=0,0,20
    r_jP = np.array([-1.026060401e-03,  2.819077873e-03])
    r_jR = np.array([1.026060401e-03, -2.819077873e-03])
    
    d = DistanceLineNode(r_iP, r_jP, r_jR=r_jR, theta_jP=45., theta_jR=45.)
    # test 1
    # d.normal = np.array([-9.998476952e-01, -1.745240630e-02])
    # test 20
    d.normal = np.array([-9.396926243e-01, -3.420201338e-01])
    # print d._distance
    
    # plt.plot([r_jP[0], r_jR[0]], [r_jP[1], r_jR[1]])
    # plt.plot(r_iP[0], r_iP[1], "o")
    
    proj = d.contact_point_on_line()
    # print "point on line =", proj
    # pprint(vars(d))
    # print '----------------------'
    print 'self._distance_sign=', d._distance_sign
    R = []
    dd = []
    # print 'r_ix=', r_ix
    for i in range(0, len(r_ix)):
        r_iPx = np.array([r_ix[i], 0.0])
        d = DistanceLineNode(r_iPx, r_jP, r_jR=r_jR, theta_jP=45., theta_jR=45.)
        dist = d._distance_sign
        R.append(r_ix[i])
        dd.append(dist)
    print 'R=', R[5]
    print 'dd=', dd[5]
    plt.plot(R, dd)
    plt.show()

    # print 'Rr=', Rr

        # plt.plot(r_ix[i],dist)


    # print 'geometry_LCS=', d.geometry_LCS()
    # print d._evaluate_distance_vector(r_iP,r_jP)
    # delta0 = -4.56499480381e-09
    # _delta = -4.5649948038044401e-09
    # _deltaF = -4.5094789549679149e-09
    # print 'delta condition=', _delta >= delta0
    # print 'delta condition=', _deltaF >= delta0

    # print 'r_iPx=',r_iPx
    # print d._evaluate_angles()
    # plt.plot(proj[0], proj[1], "s")
    # plt.show()

    
    
    
    
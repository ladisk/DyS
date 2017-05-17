__author__ = 'lskrinjar'
import time
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
from MBD_system.dr_contact_point_uP import dr_contact_point_uP


class DistanceLineNode(object):#Distance
    """
    classdocs
    """


    def __init__(self, r_iP, r_jP, normal=None, r_jR=None, tangent=None, u_iP=None, u_jP=None, normal_LCS=None, u_jR=None, tangent_LCS=None, body_id_i=None, body_id_j=None, edge_num=None, parent=None):
        """
        :param r_iP:        free a point in GCS on body i
        :param r_jP:        line point in GCS on body j
        :param normal:      normal in GCS
        :param r_jR:        line point in GCS on body j
        :param tangent:     tangent in GCS
        """
        # super(DistancePlaneNode, self).__init__(self._name, parent)
        #   parent
        self._parent = parent

        #    color
        if hasattr(self._parent, "color_GL"):
            self.color_GL = self._parent.color_GL
        else:
            self.color_GL = np.random.rand(3)

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

        #   edge number
        self.edge_num = edge_num

        #   normal, tangent in LCS
        self.tangent_LCS = tangent_LCS
        self.normal_LCS = normal_LCS
        self._contact_point_found = False

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
            #   distance
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
        fi1 = np.arccos(np.dot(self.tangent, self._distance_vector) / np.linalg.norm(self._distance_vector, ord=2))
        self._distance_vector_2 = self._evaluate_distance_vector(self.r_iP, self.r_jR)
        fi2 = np.arccos(np.dot(-self.tangent, self._distance_vector_2) / np.linalg.norm(self._distance_vector_2, ord=2))

        return fi1, fi2

    def _evaluate_distance_vector(self, r_iP, r_jP):
        """
        Evaluate distance vector from
        :return:
        """
        _distance_vector = r_iP - r_jP

        return _distance_vector

    def _evaluate_distance(self):
        """

        :return:
        """
        #   distance projection
        _distance_vector_projection = self.normal * np.dot(self._distance_vector, self.normal)

        #   distance
        self._distance = np.linalg.norm(_distance_vector_projection, ord=2)

        #    distance sign
        self._inside = False
        self._distance_sign = self._distance

        if (np.sign(np.cross(self.tangent, self._distance_vector)) >= 0) and (self.r_jR is not None):
            self.fi1, self.fi2 = self._evaluate_angles()
            # print 'self.fi1=', np.rad2deg(self.fi1), 'self.fi2=', np.rad2deg(self.fi2), self.fi1 <= np.pi/2. and self.fi2 <= np.pi/2.
            if self.fi1 <= np.pi/2. and self.fi2 <= np.pi/2.:
                self._inside = True
                self._distance_sign = -self._distance
        # print '_evaluate_distance=', self._distance_sign

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
        Function checks if distance bewteem free node and edge is inside maximum penetration depth
        :param max_penetration_depth:
        :return:
        """
        bc = self.r_iP - self.r_jP
        ba = self.r_iP - self.r_jR

        #   dot product
        dotbac = np.dot(bc, ba)

        #   theta angle
        self.theta = np.rad2deg(np.arccos(dotbac / (np.linalg.norm(bc, ord=2) * np.linalg.norm(ba, ord=2))))

        if (self.theta <= 180.) and (self._distance_sign < 0.) and ( abs(self._distance_sign) <= max_penetration_depth):
            return True

        return False

    def geometry_LCS_NOK(self, q):#NOK
        """

        :param q:
        :return:
        """
        for r, u, body_id in zip([self.r_iP, self.r_jP, self.r_jR], ["u_iP", "u_jP", "u_jR"], [self.body_id_i, self.body_id_j, self.body_id_j]):
            R = q2R_i(q, body_id)
            theta = q2theta_i(q, body_id)

            #   point coordinate in LCS
            _u = gcs2cm_lcs(r, R, theta)

            #   set object attribute
            setattr(self, u, _u)

        theta = q2theta_i(q, self.body_id_j)
        #   normal in LCS
        self.normal_LCS = gcs2cm_lcs(self.normal, np.zeros(2), theta)
        #   tangent in LCS
        self.tangent_LCS = gcs2cm_lcs(self.tangent, np.zeros(2), theta)

    def _paint_GL(self):
        """

        :return:
        """
        #   free node
        glColor3f(self.color_GL[0], self.color_GL[1], self.color_GL[2])
        glVertex3f(self.r_iP[0], self.r_iP[1], self.z_dim)

        #   edge nodes
        glVertex3f(self.r_jP[0], self.r_jP[1], self.z_dim)
        glVertex3f(self.r_jR[0], self.r_jR[1], self.z_dim)


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
    r_iP = np.array([.5, .6]) #OK
    r_iP = np.array([0.011439999485,  0.006928203106]) #NOK
    #   line node
    r_jP = np.array([0.004, -0.006928203106])
    r_jR = np.array([8.000000000000e-03,  -1.959434964136e-18])

    #   free node
    # r_iP = np.array([-0.007592, -0.011551])
    # # r_iP = np.array([.4, .6])
    # #   line node
    # r_jR = np.array([-0.007, -0.01])
    # r_jP = np.array([-0.009, -0.015])
    
    d = DistanceLineNode(r_iP, r_jP, r_jR=r_jR)
    
    print d._distance
    
    plt.plot([r_jP[0], r_jR[0]], [r_jP[1], r_jR[1]])
    plt.plot(r_iP[0], r_iP[1], "o")
    
    proj = d.contact_point_on_line()
    print "point on line =", proj
    pprint(vars(d))
    plt.plot(proj[0], proj[1], "s")
    plt.show()

    
    
    
    
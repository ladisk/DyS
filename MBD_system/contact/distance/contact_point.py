"""

created by: lskrinjar
date of creation: 21/07/2016
time of creation: 21:55
"""
import os
import inspect
import time
from pprint import pprint
import numpy as np
from matplotlib import pyplot as plt


from MBD_system.contact.distance.distance_line_node import DistanceLineNode
from MBD_system.transform_cs import gcs2cm_lcs, cm_lcs2gcs, uP_lcs2gcs
from MBD_system.q2R_i import q2R_i
from MBD_system.q2theta_i import q2theta_i
from MBD_system.q2dR_i import q2dR_i
from MBD_system.q2dtheta_i import q2dtheta_i
from MBD_system.n2t import n2t
from MBD_system.dr_contact_point_uP import dr_contact_point_uP
from MBD_system.force.force import Force


class ContactPoint(DistanceLineNode):
    """
    classdocs
    """

    def __init__(self, dist_obj=None, parent=None):
        """
        :param r_iP:        free a point in GCS on body i
        :param r_jP:        line point in GCS on body j
        :param normal:      normal in GCS
        :param r_jR:        line point in GCS on body j
        :param tangent:     tangent in GCS
        """
        # stack = inspect.stack()
        # the_class = stack[1][0].f_locals["self"].__class__
        # the_method = stack[1][0].f_code.co_name
        # print("  I was called by {}.{}()".format(str(the_class), the_method))

        super(ContactPoint, self).__init__(np.zeros(2), np.zeros(2), parent=parent)
        #   assign all distance object attributes to this object attributes
        if hasattr(dist_obj, "__dict__"):
            self.__dict__ = dist_obj.__dict__

            #   overlap pair
            self.overlap_pair = self._parent
            # print "self.overlap_pair@__init__ =", self.overlap_pair

            #   parent
            self._parent = parent

        #   status
        self._contact_point_found = False
        self.active = True

        #   coordinates in GCS
        self.r_P_GCS_list = [np.zeros(2, dtype="float32"), np.zeros(2, dtype="float32")]
        self.u_P_LCS_list = [np.zeros(2, dtype="float32"), np.zeros(2, dtype="float32")]

        #   normals and tangents in GCS
        self._n_GCS_list = [np.zeros(2, dtype="float32"), np.zeros(2, dtype="float32")]
        self._t_GCS_list = [np.zeros(2, dtype="float32"), np.zeros(2, dtype="float32")]

        #   contact point velocity
        self._dq_n = 0.
        self._dq_t = 0.
        self._dq0_n = 0.
        self._dq0_t = 0.

        #   contact force
        self.Fn = 0.
        self.Ft = 0.
        self.F = np.zeros(2)
        self._Fn_list = self._create_forces(direction="normal")
        self._Ft_list = self._create_forces(direction="tangent")

    def _create_forces(self, direction=""):
        """
        Method creates forces objects for each body in contact
        :return:
        """
        F_list = []
        if hasattr(self._parent, "body_id_list"):
            #   create force object
            if direction == "normal":
                dir = "n"
            elif direction == "tangent":
                dir = "t"
            else:
                dir = direction

            # for body_id, visibility in zip(self._parent.body_id_list, self._parent._visible_F_list):
            for body_id, visibility in zip(self.body_id_list, self._parent._visible_F_list): #  corrected from self._parent.body_id_list, self._parent._visible_F_list to self.body_id_list, self._parent._visible_F_list
                #   create force object
                F = Force(body_id, force_name=self._parent._name + "_F" + dir + "_on_body_" + str(body_id), visible=visibility, parent=self._parent)
                F_list.append(F)

                #   create force visualization data
                ren = self._parent._parent._parent._parent.dys.simulation_control_widget.vtkWidget.renderer
                F.set_vtk_data(renderer=ren)

                # if F.vtk_actor is not None:
                #     self._parent._parent._parent._parent.dys.simulation_control_widget.vtkWidget.renderer.AddActor(F.vtk_actor)
                # else:
                #     print Warning, "Force visualizations not created!"

                #   append force object to list of forces of MBD system object
                self._parent._parent._parent.forces.append(F)

                # #   append force object to list of forces of body object
                # self._parent._parent._parent.bodies[body_id].forces.append(F)

        return F_list

    def deactivate_forces(self):
        """

        :return:
        """
        print "deactivate_forces()@",__name__
        #   join both lists to one
        F_list = self._Fn_list + self._Ft_list

        for i, F in enumerate(F_list):
            print "i =", i, F._name
            F.deactivate()

            # F._clear_vtk_data()

        self._Fn_list = []
        self._Ft_list = []
        # for i, F in enumerate(F_list):
        #     self._parent._parent._parent.forces.remove(F)
        print "N@deactivate_forces() =", len(self._parent._parent._parent.forces)
        # self.active = False

    def remove_forces(self):
        """

        :return:
        """
        F_list = self._Fn_list + self._Ft_list

        for F in F_list:
            F.clear_vtk_data()

        # for F_list in F_nt_list:
        #     self._remove_forces(F_list)
        #
        # for F_list in F_nt_list:
        #     self._delete_forces(F_list)

    # def _delete_forces(self, F_list):
    #     """
    #     Delete force objects
    #     :return:
    #     """
    #     print "_delete_forces()"
    #     for F in F_list:
    #         del F

    # def _remove_forces(self, F_list):
    #     """
    #     Remove force objects from lists
    #     :return:
    #     """
    #     print "_remove_forces()"
    #     for F, body_id in zip(F_list, self._parent.body_id_list):
    #         if F in self._parent._parent._parent.forces:
    #             self._parent._parent._parent.forces.remove(F)
    #
    #         if F in self._parent._parent._parent.forces.bodies[body_id].forces:
    #             self._parent._parent._parent.forces.bodies[body_id].forces.remove(F)

    def _clear_vtk_data(self):
        """

        :param F_list:
        :return:
        """
        F_list = self._Fn_list + self._Ft_list
        for F in F_list:
            if F.vtk_actor is not None and F.renderer is not None:
                F._clear_vtk_data()

    # def _remove_forces_TODO(self, F_list):
    #     """
    #
    #     :param F_list:
    #     :return:
    #     """
    #     print "_remove_forces()@",__name__
    #     if hasattr(self._parent, "body_id_list"):
    #         for F, body_id in zip(F_list, self.body_id_list):
    #             #   remove from renderer
    #             if F.vtk_actor is not None:
    #
    #                 self._parent._parent._parent._parent.dys.simulation_control_widget.vtkWidget.renderer.RemoveViewProp(F.vtk_actor)
    #                 self._parent._parent._parent._parent.dys.simulation_control_widget.vtkWidget.renderer.RemoveActor(F.vtk_actor)
    #
    #                 F.force.Delete()
    #                 F.vtk_mapper.Delete()
    #
    #             else:
    #                 print Warning, "Force visualizations not removed!"
    #
    #             #   append force object to list of forces of MBD system object
    #             self._parent._parent._parent.forces.remove(F)
    #
    #             #   append force object to list of forces of body object
    #             self._parent._parent._parent.bodies[body_id].forces.remove(F)
    #
    #             del F

    def evaluate_contact_points(self, q):
        """
        Method evaluated actual contact point:
        on body i is actually already evaluated as free node
        on body j is a projection on edge
        First in GCS and than in LCS and attributes r_P_GCS_list are u_P_LCS_list are evaluated
        :return:
        """
        #   Checked and working correctly - OK!
        #   in GCS
        self.r_P_GCS_list[0] = self.r_iP
        self.r_P_GCS_list[1] = self.contact_point_on_line()

        self._contact_geometry_LCS(q)

    def _contact_geometry_LCS(self, q):
        """

        :param q:
        :return:
        """
        #   in LCS
        # Checked, the method executed correctly - OK!

        for i, (body_id, rP) in enumerate(zip(self.body_id_list, self.r_P_GCS_list)):
            R = q2R_i(q, body_id)
            theta = q2theta_i(q, body_id)
            uP = gcs2cm_lcs(rP, R, theta)

            self.u_P_LCS_list[i] = uP

    def update_contact_points_GCS(self, q):
        """
        Function updates coordinates of contact points (geometry) in GCS based on LCS coordinates and vector q
        :param q: vector of MBD system coordinates
        :return:
        """
        # print "update_contact_points_GCS@", __file__
        for i, (body_id, uP) in enumerate(zip(self.body_id_list, self.u_P_LCS_list)):
            R = q2R_i(q, body_id)
            theta = q2theta_i(q, body_id)

            #   point coordinate in LCS
            rP = cm_lcs2gcs(uP, R, theta)

            self.r_P_GCS_list[i] = rP

        #   distance and delta are equal and multiplied with (-1)
        self._distance_vector = self._evaluate_distance_vector(self.r_P_GCS_list[0], self.r_P_GCS_list[1])

        # _distance = self._distance_sign = - np.linalg.norm(self._distance_vector, ord=2)
        _distance = self._distance_sign = - np.linalg.norm((self.normal * np.dot(self._distance_vector, self.normal)), ord=2)
        return _distance, self._distance_sign, self.normal, self.tangent

    def contact_geometry(self, q):
        """
        Normals and tangens in LCS and GCS of each body
        :return:
        """
        #   contact normals and tangents in LCS
        theta = q2theta_i(q, self.body_id_j)
        #   normal in LCS
        n_LCS = gcs2cm_lcs(self.normal, np.zeros(2), theta)
        self._n_LCS_list = [n_LCS, -n_LCS]
        #   tangent in LCS
        t_LCS = gcs2cm_lcs(self.tangent, np.zeros(2), theta)
        self._t_LCS_list = [t_LCS, -t_LCS]

        # for i, (body_id, n_LCS_i, t_LCS_i) in enumerate(zip(self.body_id_list, self._n_LCS_list, self._t_LCS_list)):
        #     theta = q2theta_i(q, body_id)
        #     self._n_GCS_list[i] = uP_lcs2gcs(n_LCS_i, theta)
        #     self._t_GCS_list[i] = uP_lcs2gcs(t_LCS_i, theta)

        for i, (body_id, sign) in enumerate(zip(self.body_id_list, [+1., -1.])):
            theta = q2theta_i(q, body_id)

            self._n_GCS_list[i] = self.normal * sign
            self._t_GCS_list[i] = self.tangent * sign

    def _contact_velocity(self, q):
        """

        Args:
            q:

        Returns:

        """
        dr_P = []
        for body_id, u_P in zip(self.body_id_list, self.u_P_LCS_list):
            #   body velocity, R, theta
            dR = q2dR_i(q, body_id)
            theta = q2theta_i(q, body_id)

            #    dtheta - omega
            dtheta = q2dtheta_i(q, body_id)

            #    point velocity
            dr_P_body = dr_contact_point_uP(dR, theta, dtheta, u_P)
            # print 'dr_P_body=', dr_P_body, 'body_id=', body_id

            #    add to list
            dr_P.append(dr_P_body)

        _dq = dr_P[1] - dr_P[0]

        self._dq_n = np.dot(_dq, self.normal)
        # print 'dq_n=', self._dq_n, 'normal=', self.normal, _dq
        self._dq_t = np.dot(_dq, self.tangent)

        return self._dq_n, self._dq_t

    def set_dq0(self, dq0_n=None, dq0_t=None):
        """
        Function saves the initial normal and tangential contact velocity at impact
        :param dq0_n    a normal contact velocity (scalar, float value)
        :param dq0_t    a tangential contact velocity (scalar, float value)
        """
        if dq0_n is None:
            self._dq0_n = self._dq_n
        else:
            self._dq0_n = dq0_n

        if dq0_t is None:
            self._dq0_t = self._dq_t
        else:
            self._dq0_t = dq0_t

    def evaluate_F(self, Fn, Ft):
        """

        Returns:

        """
        self.Fn = Fn
        self.Ft = Ft
        self.F = self.Fn * self.normal + self.Ft * self.tangent

    def __del__(self):
        """

        :return:
        """
        # self._remove_forces(self._Fn_list)
        # self._remove_forces(self._Ft_list)


if __name__ == "__main__":
    contact_point = ContactPoint()

    # rjP = np.array([0, 3E-3])
    # rjR = np.array([0, -3E-3])
    # n = np.array([-1, 0])

    riP = np.array([-8.121906820e-03,  1.128993071e-03])

    rjP = np.array([-8.086623670e-03,  0.000000000e+00])
    rjR = np.array([-8.186623671e-03,  3.200000048e-03])
    n = np.array([9.995120761e-01,  3.123475238e-02])

    contact_point.r_iP = riP
    contact_point.r_jP = rjP
    contact_point.r_jR = rjR
    contact_point.theta_jP = 45.
    contact_point.theta_jR = 45.

    pprint(vars(contact_point))

    # riP = np.array([0, 2E-4])
    contact_point.normal = n

    # riPx = np.linspace(0, 1E-4, 100)
    # delta = np.zeros_like(riPx)

    # for i in xrange(0, len(riPx)):
    #     riP[0] = riPx[i]
    #     contact_point.r_iP = riP
    #
    #     contact_point._distance_vector = contact_point._evaluate_distance_vector(contact_point.r_iP, contact_point.r_jP)
    #     contact_point._evaluate_distance()
    #
    #     delta[i] = contact_point._distance_sign
    #     print riPx[i], delta[i]


    contact_point.plot()
    plt.show()




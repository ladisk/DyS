"""

created by: lskrinjar
date of creation: 25/07/2016
time of creation: 21:40
"""
from pprint import pprint
import itertools
import numpy as np
import sympy as sp
import vtk
from matplotlib import pyplot as plt

from MBD_system.ancf.mesh.cross_section import cross_section

np.set_printoptions(precision=4, threshold=None, edgeitems=100, linewidth=1000, suppress=False, nanstr=None, infstr=None)


from MBD_system.body.body import Body
from MBD_system.ancf.mesh.finite_element.beam_2d import Beam2D
from MBD_system.q2R_i import q2R_i
from MBD_system.gaussian_quadrature_on_arbitrary_intervals import change_interval
from MBD_system.ancf.mesh.cross_section.cross_section import CrossSection
from MBD_system.A import A_matrix


class Beam2DEulerBernoulli(Beam2D):
    """
    Two dimensional Euler-Bernoulli beam element
    Computational Continuum mechanics 2008, A. A. Shabana

    Each node k (k=i, j) has four coordinates; two translations rjk, and two gradient
    coordinates rjk_x1. Therefore, the vector of nodal coordinates has eight elements.
    """

    def __init__(self, id, node_id_i, node_id_j, mesh=None, nodes=None, parent=None):
        """

        """
        super(Beam2DEulerBernoulli, self).__init__(id, node_id_i, node_id_j, mesh=mesh, nodes=nodes, parent=parent)

        #   body coordinates
        self.q_i_size = Beam2DEulerBernoulli.get_e_size()

        #   element type
        self.finite_element_type = Beam2DEulerBernoulli.get_element_type()

        #   vector of nodal coordinates
        #   vector of time-dependent coefficients or
        #   coordinates, which consist of absolute position and slope coordinates (see Ch. 1 and 2)
        self.e_n = Beam2DEulerBernoulli.get_e_size()
        #   number of absolute coordinates per element node
        self.node_e_n = self.e_n/2

        self.e0 = self.evaluate_e0()
        self.e = None
        #   vector of velocities
        self.de = np.zeros(self.e_n)

        #   order of shape function polynomials
        self.N = 2

        #   size of mass matrix
        self.M_size = self.e_n

        #   predefine zero matrix size for each node
        self._T_ijk_cols = 6#or 4 - todo

        #   gradient deficient properties
        self.g_1 = np.array([1, 0], dtype="float")
        self.A = A_matrix(self.theta[2])

        #   constraint matrix
        self.C_q = None

        #   predefined matrices attributes
        self.K_t = None
        self.K_l = None

        if self.mesh is not None:
            #   cross-section
            if self.mesh.cross_section is not None:
                self.Izz = self.mesh.cross_section.Izz
                self.area = self.mesh.cross_section.area

            if hasattr(self.mesh, "_parent"):
                if hasattr(self.mesh._parent, "module_of_elasticity"):
                    self.module_of_elasticity = self.mesh._parent.module_of_elasticity

                elif hasattr(self.mesh, "module_of_elasticity"):
                    self.module_of_elasticity = self.mesh.module_of_elasticity

                else:
                    pass

                if hasattr(self.mesh._parent, "density"):
                    self.density = self.mesh._parent.density

    @staticmethod
    def get_e_size():
        """

        :return:
        """
        return 8

    @staticmethod
    def get_element_type():
        """

        :return:
        """
        return "beam 2D euler-bernoulli"

    def evaluate_e0(self):
        """
        Vector of absolute position and slope coordinates
        :return:
        """
        e = np.zeros(self.e_n)

        #   coordinates
        if self.mesh is not None:
            e[0:2] = self.mesh.nodes0[self.node_id_i][0:2]
            e[4:6] = self.mesh.nodes0[self.node_id_j][0:2]

        #   gradient coordinates
        de = np.array([np.cos(self.theta[2]), np.sin(self.theta[2])])
        e[2:4] = de
        e[6:8] = de

        return e

    def evaluate_e(self):
        """
        Vector of absolute position and slope coordinates
        :return:
        """
        e = np.zeros(self.e_n)

        #   coordinates
        e[0:2] = self.mesh.nodes[self.node_id_i][0:2]
        e[4:6] = self.mesh.nodes[self.node_id_j][0:2]

        #   beam geometry
        #   length of the beam
        self.L_vector = e[4:6] - e[0:2]
        self.L = np.linalg.norm(self.L_vector)

        #   theta
        self.theta[2] = np.arctan2(self.L_vector[1], self.L_vector[0])

        #   gradient vectors
        de = np.array([np.cos(self.theta[2]), np.sin(self.theta[2])])
        e[2:4] = de
        e[6:8] = de

        return e

    def evaluate_q(self):
        """
        Function returns vector of absolute nodal coordinates of finite element
        :return:
        """
        return self.e

    def evaluate_dq(self):
        """

        :return:
        """
        return self.de

    def evaluate_r(self, ksi, e_i=None):
        """
        Function evaluates coordinates in GCS based on LCS data and shape function of finite element
        :param uPj:
        :return:
        """
        S = self._evaluate_S(ksi)

        if e_i is None:
            if self.e is None:
                self.e = self.evaluate_e()

            e_i = self.e

        r = np.dot(S, e_i)
        return r

    def _evaluate_S(self, ksi):
        """

        :return:
        """
        #   shape functions
        S1 = self._S1(ksi)
        S2 = self._S2(ksi)
        S3 = self._S3(ksi)
        S4 = self._S4(ksi)

        S = np.hstack(([S1 * np.eye(2), S2 * np.eye(2), S3 * np.eye(2), S4 * np.eye(2)]))
        return S

    def _S1(self, ksi):
        """

        :param ksi:
        :return:
        """
        S1 = 1. - 3. * ksi ** 2 + 2. * ksi ** 3
        return S1

    def _S2(self, ksi):
        """

        :param ksi:
        :return:
        """
        S2 = self.L * (ksi - 2. * ksi**2 + ksi**3)
        return S2

    def _S3(self, ksi):
        """

        :param ksi:
        :return:
        """
        S3 = 3. * ksi**2 - 2. * ksi**3
        return S3

    def _S4(self, ksi):
        """

        :param ksi:
        :return:
        """
        S4 = self.L * (ksi ** 3 - ksi ** 2)
        return S4

    def _evaluate_S_symbolic(self, x=None):
        """

        :param x:
        :return:
        """
        ksi, L = sp.symbols('ksi L')

        S1 = sp.integrate(1 - 3. * ksi**2 + 2. * ksi**3, (ksi, 0, 1)).subs({L: self.L})
        S2 = sp.integrate(L * (ksi - 2. * ksi**2 + ksi**3), (ksi, 0, 1)).subs({L: self.L})
        S3 = sp.integrate(3. * ksi**2 - 2. * ksi**3, (ksi, 0, 1)).subs({L: self.L})
        S4 = sp.integrate(L * (-ksi**2 + ksi**3), (ksi, 0, 1)).subs({L: self.L})

        S = np.hstack(([S1 * np.eye(2), S2 * np.eye(2), S3 * np.eye(2), S4 * np.eye(2)]))
        return S

    def _evaluate_dS_dksi(self, ksi):
        """

        :return:
        """
        #   shape functions
        dS1 = -6. * ksi + 6. * ksi**2
        dS2 = self.L * (1. - 4. * ksi + 3. * ksi**2)
        dS3 = 6. * ksi - 6. * ksi**2
        dS4 = self.L * (-2. * ksi + 3. * ksi**2)

        dS = np.hstack(([dS1 * np.eye(2), dS2 * np.eye(2), dS3 * np.eye(2), dS4 * np.eye(2)]))
        return dS

    def _evaluate_ddS_ddksi(self, ksi):
        """

        :return:
        """
        #   shape functions
        ddS1 = self._ddS1(ksi)
        ddS2 = self._ddS2(ksi)
        ddS3 = self._ddS3(ksi)
        ddS4 = self._ddS4(ksi)

        ddS = np.hstack(([ddS1 * np.eye(2), ddS2 * np.eye(2), ddS3 * np.eye(2), ddS4 * np.eye(2)]))
        return ddS

    def _ddS1(self, ksi):
        ddS1 = -6. + 12. * ksi
        return ddS1

    def _ddS2(self, ksi):
        ddS2 = (-4 + 6. * ksi) * self.L
        return ddS2

    def _ddS3(self, ksi):
        ddS3 = 6. - 12. * ksi
        return ddS3

    def _ddS4(self, ksi):
        ddS2 = (-2. + 6. * ksi) * self.L
        return ddS2

    def evaluate_M(self):
        """

        :return:
        """
        M = np.zeros([self.e_n, self.e_n])

        #   get coordinates and weight of legendre polynomials for gauss quadrature
        x, w = np.polynomial.legendre.leggauss(self.N * 2)

        for x_i, w_i in zip(x, w):
            x_i_interval, ab_i_interval = change_interval(x_i, self.a, self.b)
            M_i = w_i * np.dot(self._evaluate_S(x_i_interval).T, self._evaluate_S(x_i_interval)) * ab_i_interval

            M += M_i

        if self.mass == 0.:
            raise Warning, "Element attribute mass is zero. Element mass matrix is zero matrix."

        return M * self.mass

    def evaluate_M_testing(self):
        """

        :return:
        """
        M = np.zeros([self.e_n, self.e_n])

        #   get coordinates and weight of legendre polynomials for gauss quadrature
        x, w = np.polynomial.legendre.leggauss(self.N * 2)

        for x_i, w_i in zip(x, w):
            x_i_interval, ab_i_interval = change_interval(x_i, self.a, self.b)
            M_i = w_i * np.dot(self._evaluate_S(x_i_interval).T, self._evaluate_S(x_i_interval)) * ab_i_interval

            M += M_i

        if self.mass == 0.:
            raise Warning, "Element attribute mass is zero. Element mass matrix is zero matrix."

        else:
            M = M * self.mass

        return M

    def evaluate_C(self):
        """
        Function evaluates damping matrix of a beam element
        :return:
        """
        C = np.zeros([self.e_n, self.e_n])

        #   get coordinates and weight of legendre polynoms for gauss quadrature
        x, w = np.polynomial.legendre.leggauss(self.N * 2)

        self.C_j = ((self.b - self.a) / 2.) * C
        return self.C_j

    def evaluate_Q_k(self, e):
        """

        :param q:
        :return:
        """
        Q_j = np.dot(self.K_j, e)
        return Q_j

    def evaluate_Q_s(self, e):
        """

        :param q:
        :return:
        """
        K = self.evaluate_K(e)
        Q_s = np.dot(K, e)
        return Q_s

    def evaluate_Q_e_M(self, e_i, M, ksi):
        """

        :param e_i:
        :param M:   moment (of type float)
        :param ksi: ksi
        :return:
        """
        # print "evaluate_Q_e_M()",__name__
        # print "ksi =", ksi
        dS_dksi = self._evaluate_dS_dksi(ksi)
        # print "np.dot(dS_dksi, e_i) =", np.dot(dS_dksi, e_i)
        d = np.linalg.norm(np.dot(dS_dksi, e_i), ord=1)
        # print "d1 =", d1
        # d2 = e_i[2]**2 + e_i[3]**2
        # print "d2 =", d2

        Q_e_M = (M / (d * self.L**2)) * ((np.dot(dS_dksi[0], e_i) * dS_dksi[1]) - (np.dot(dS_dksi[1], e_i) * dS_dksi[0]))
        # print "Q_e_M1 =", Q_e_M

        # Q_e_M = np.zeros_like(e_i)
        # if ksi == 0:
        #     Q_e_M[2] = -(M * e_i[3]) / d2
        #     Q_e_M[3] = (M * e_i[2]) / d2
        #
        # if ksi == 1:
        #     Q_e_M[6] = -(M * e_i[7]) / d2
        #     Q_e_M[7] = (M * e_i[6]) / d2
        # print "Q_e_M2 =", Q_e_M

        return Q_e_M

    def evaluate_K(self, e_i):
        """
        See ref: Development of simple models for the elastic forces in the absolute nodal co-ordinate formulation (doi:10.1006/jsvi.1999.2935)
        :param e:
        :return:
        """
        self.K = self._evaluate_K_t() + self._evaluate_K_l(e_i)
        return self.K

    def _construct_constant_matrix_K_t(self):
        """

        :return:
        """
        #   upper diagonal matrix
        K_u = np.array([[12., 0., 6.*self.L, 0., -12., 0., 6.*self.L, 0.],
                        [0., 12., 0., 6.*self.L, 0., -12., 0., 6.*self.L],
                        [0., 0., 4.*self.L**2, 0., -6.*self.L, 0., 2.*self.L**2, 0.],
                        [0., 0., 0., 4.*self.L**2, 0., -6.*self.L, 0., 2.*self.L**2],
                        [0., 0., 0., 0., 12., 0., 6*self.L, 0.],
                        [0., 0., 0., 0., 0., 12., 0., -6.*self.L],
                        [0., 0., 0., 0., 0., 0., 4.*self.L**2, 0.],
                        [0., 0., 0., 0., 0., 0., 0., 4.*self.L**2]])

        K = (K_u + K_u.T - np.diag(K_u.diagonal())) * (self.module_of_elasticity * self.mesh.cross_section.Izz / self.L**3)

        return K

    def _evaluate_K_t(self):
        """
        Function evaluates stiffness matrix of a beam element
        :return:
        """
        if self.K_t is None:
            K = np.zeros([self.e_n, self.e_n])

            #   get coordinates and weight of legendre polynomials for gauss quadrature
            x, w = np.polynomial.legendre.leggauss(self.N * 2)

            for x_i, w_i in zip(x, w):
                x_i_interval = ((((self.b - self.a) * x_i) + (self.b + self.a)) / 2.)
                K_i = w_i * np.dot(self._evaluate_ddS_ddksi(ksi=x_i_interval).T, self._evaluate_ddS_ddksi(ksi=x_i_interval))

                K += K_i

            #   stiffness matrix of i-th element
            self.K_t = ((self.b - self.a) / 2.) * (self.module_of_elasticity * self.Izz / (self.L**3)) * K

        return self.K_t

    def _evaluate_K_l(self, e_i):
        """
        See ref: development of simple models for the elastic forces in the absolute nodal co-ordinate formulation (doi:10.1006/jsvi.1999.2935)
        :param e:
        :return:
        """
        #   distance between the nodes
        d = self.evaluate_L(e_i)

        #   average longitudinal strain along the element
        epsilon_l = (d - self.L) / self.L

        #   predefine size
        K = np.zeros([self.e_n, self.e_n])

        #   get coordinates and weight of legendre polynomials for gauss quadrature
        x, w = np.polynomial.legendre.leggauss(self.N * 2)

        for x_i, w_i in zip(x, w):
            x_i_interval = ((((self.b - self.a) * x_i) + (self.b + self.a)) / 2.)
            K_i = w_i * np.dot(self._evaluate_dS_dksi(ksi=x_i_interval).T, self._evaluate_dS_dksi(ksi=x_i_interval))

            K += K_i

        #   stiffness matrix of i-th element
        K = (self.module_of_elasticity * self.area / self.L) * epsilon_l * K

        return K

    def evaluate_B(self):
        """
        Only nodal coordinates vector of first node of every element is added to the vector of nodal coordinates of the mesh.
        And from the last element in a mesh also the absolute nodal coordinates of the last node is also added to the nodal coordinates of the mesh.

        Boolean matrix of the finite element always has a number of rows
        equal to the number of the finite element nodal coordinates and a number of
        columns equal to the number of the body nodal coordinates.
        :return:
        """
        if self.B is None:
            #   predefine zero matrix
            e_n = sum(self.mesh.node_dim[self.node_id_i:self.node_id_j+1])
            self.B = np.zeros([e_n, self.mesh.n_NC])

            if self.node_id_list[0] == 0:
                _i = 0
            else:
                _i = sum(self.mesh.node_dim[0:self.node_id_list[0]])

            # print "first node id in element =", self.node_id_list[0]
            # print "index start =", self.mesh.node_dim[self.node_id_list[0] - 1]
            # print "index end =", self.mesh.node_dim[self.node_id_list[0]]+e_n
            self.B[:, _i:_i+e_n] = np.eye(e_n)
            # print "self.mesh.node_dim =", self.mesh.node_dim
            # for i, element in enumerate(self.mesh.elements):
            #     if element is not self:
            #         if self.node_id_i in element.node_id_list:
            #             print "i@evaluate_B() =", i
            #             self.B[:, (e_n/2) * (self.node_id_i + 1):(e_n/2) * (self.node_id_i + 1) + e_n] = np.eye(e_n)
            #
            #         if self.node_id_j in element.node_id_list:
            #             self.B[:, (e_n/2) * (self.node_id_j - 1):(e_n/2) * (self.node_id_j - 1) + e_n] = np.eye(e_n)
            #
            #     if len(self.mesh.elements) == 1:
            #         self.B = np.eye(self.e_n)

        return self.B

    def evaluate_T(self):
        """

        :return:
        """
        if self.T is None:
            self.T = np.zeros([self.e_n, sum([self.mesh.node_dim[self.node_id_i], self.mesh.node_dim[self.node_id_j]])])

            for i, (node_id, ksi) in enumerate(zip(self.node_id_list, self._ksi_list)):
                #   transformation matrix of a node on element
                T_ijk = self._evaluate_T_ijk(ksi, node_id)
                # print T_ijk
                if node_id == self.node_id_i:
                    self.T[self.node_e_n * i:self.node_e_n * i + self.node_e_n, 0:self.mesh.node_dim[self.node_id_i]] = T_ijk
                if node_id == self.node_id_j:
                    self.T[self.node_e_n * i:self.node_e_n * i + self.node_e_n, self.mesh.node_dim[self.node_id_i]:self.mesh.node_dim[self.node_id_i]+self.mesh.node_dim[self.node_id_j]] = T_ijk

        return self.T

    def _evaluate_T_ijk(self, ksi, node_id):
        """
        Transformation matrix of a node
        :param ksi: parameter ksi
        :return:
        """
        if self.mesh.node_dim[node_id] == 6:
            #   predefine matrix size
            T = np.zeros([self.node_e_n, self._T_ijk_cols])
            T[0:2, 0:2] = np.eye(2, dtype="float")

            #   value of derivative of shape function matrix against parameter ksi
            # dSksi = self._evaluate_dS_dksi(ksi)

            #   elements of matrix T_ijk
            # j_11 = np.dot(dSksi[0, :], self.e0)
            # j_21 = np.dot(dSksi[1, :], self.e0)
            # print "j_11 =", j_11
            # print "j_21 =", j_21

            # T[2:4, 2:4] = j_11 * np.eye(2, dtype="float")
            # T[2:4, 4:6] = j_21 * np.eye(2, dtype="float")
            A11 = self.A[0, 0]
            # # print "A11 =", A11
            A21 = self.A[1, 0]
            # # print "A21 =", A21
            T[2:4, 2:4] = A11 * np.eye(2, dtype="float")
            T[2:4, 4:6] = A21 * np.eye(2, dtype="float")

            #   new formulation
            # T[2:4, 2:4] = np.dot(self.A, self.g_1)
            # T[2:4, 4:6] = np.dot(self.A, self.g_1)

        elif self.mesh.node_dim[node_id] == 4:
            T = np.eye(self.node_e_n)

        else:
            T = None
            print Warning, "Transformation matrix T_ijk not constructed!"

        return T

    def evaluate_TB(self):
        """

        :return:
        """

    def evaluate_C_q_fixed(self, node_id):
        """

        :param q:
        :return:
        """
        #   node vector of position coordinates
        ksi = self.node_id_list.index(node_id)

        #   predefine matrix
        matrix = np.zeros([4, self.e_n])

        #   constrain on absolute nodal coordinates e1, e2 (x, y)
        matrix[0:2, :] = self._evaluate_S(ksi)

        #   constrain on absolute nodal coordinates e3, e4 (gradients)
        if self.node_id_list.index(node_id) == 0:
            matrix[2, 2] = 1.
            matrix[3, 3] = 1.
        if self.node_id_list.index(node_id) == 1:
            matrix[2, 6] = 1.
            matrix[3, 7] = 1.

        return matrix

    def evaluate_C_q_hinged(self, node_id):
        """

        :param q:
        :return:
        """
        #   node vector of position coordinates
        x = self.mesh.nodes[node_id]

        #   constrain on absolute nodal coordinates e1, e2 (x, y)
        matrix = np.zeros([2, self.e_n])

        #   constrain on absolute nodal coordinates e3, e4 (gradients)
        if self.node_id_list.index(node_id) == 1:
            matrix[0, 0] = 1.
            matrix[1, 1] = 1.
        if self.node_id_list.index(node_id) == 0:
            matrix[0, 4] = 1.
            matrix[1, 5] = 1.
        return matrix

    def evaluate_Q_d_fixed(self, node_id):
        """

        :param node_id:
        :return:
        """
        Q_d = 0
        return Q_d


if __name__ == "__main__":
    #   example data from: A geometrically exact beam element based on the absolute nodal coordinate formulation
    g = np.array([0., -9.83])

    m = 1.
    nodes = np.array([[0., 0.],
                      [1., 0.]])
    element = Beam2DEulerBernoulli(0, 0, 1, mesh=None, nodes=nodes)
    element.h = 1.
    element.w = 1.
    element.m = m
    cross_section = CrossSection(type="rectangle")
    cross_section.B = element.w
    cross_section.H = element.h
    cross_section.evaluate()

    element.Izz = 1.
    element.density = 1.
    element.module_of_elasticity = 1.
    print "K ="
    K = element.evaluate_K()
    print K
    print K == K.T
    # print "L =", element.L
    # print "K ="
    # print element.evaluate_K()
    # print "Q_g"
    # print "numerical =", m * element.evaluate_Q_g(g)
    # print "analytical =", m * np.dot(g, element._evaluate_S_symbolic())
    # print element._evaluate_dSdx1(np.array([0, 0]))
    # print element.evaluate_geometry(element.evaluate_e())
    element.e0 = element.e = np.array([1.5, 1.5, -3, -2, -1, 1, -2, 3], dtype="float")
    element.L = 5.
    fig = plt.figure(figsize=(6, 5),
                     dpi=100,
                     facecolor='w',
                     edgecolor='k')

    ax = plt.subplot(111, aspect="equal")
    ax.ticklabel_format(style='sci', axis='both')

    element.plot(ax)
    plt.show()
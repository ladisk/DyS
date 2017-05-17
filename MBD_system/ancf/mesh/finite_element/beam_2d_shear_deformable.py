"""

created by: lskrinjar
date of creation: 25/07/2016
time of creation: 21:40
"""
from pprint import pprint
import itertools
import numpy as np
import vtk


from MBD_system.ancf.mesh.finite_element.beam_2d import Beam2D


class Beam2DShearDeformable(Beam2D):
    """
    Two dimensional Euler-Bernoulli beam element
    Computational Continuum mechanics 2008, A. A. Shabana

    Each node k (k=i, j) has six degrees of freedome e.g., coordinates; two translations rjk, and four gradient
    coordinates defined by two vectors rjk_x1 and rjk_x2. Therefore, the vector of nodal coordinates has eight elements.
    """

    def __init__(self,  id, node_i, node_j, parent=None):
        super(Beam2DShearDeformable, self).__init__(id, node_i, node_j, parent=parent)

        #   body coordinates
        self.q_i_size = Beam2DShearDeformable.get_e_size()

        #   element type
        self.finite_element_type = Beam2DShearDeformable.get_element_type()

        #   vector of nodal coordinates
        self.e = Beam2DShearDeformable.get_e_size()

        #   order of shape function polynomials
        self.N = 4

        #   size of mass matrix
        self.M_size = self.e_n

    @staticmethod
    def get_e_size():
        """

        :return:
        """
        return 12

    @staticmethod
    def get_element_type():
        """

        :return:
        """
        return "beam 2D shear deformable"

    def evaluate_M(self, q):
        """

        :return:
        """
        M = np.zeros([self.e_n, self.e_n])

        #   get coordinates and weight of legendre polynoms for gauss quadrature
        x, w = np.polynomial.legendre.leggauss(self.N)

        for x_i, w_i in zip(x, w):
            x_i_interval = ((self.b - self.a) * x_i + (self.b + self.a)) / 2.
            M_i = w_i * np.dot(self._shapeFunction(x_i_interval).T, self._shapeFunction(x_i_interval)) * (self.b - self.a) / 2.

            M += M_i

        return self.Izz * M

    def evaluate_M_size(self):
        """
        Function evaluates and returns size of mass matrix of body
        :return:
        """
        return self.M_size

    def _evaluate_S(self, x):
        """

        :return:
        """
        #   normed length of the beam element
        ksi = x1 / self.L
        ni = x2 / self.L

        #   shape functions
        S1 = self._S1(ksi, ni)
        S2 = self._S2(ksi, ni)
        S3 = self._S3(ksi, ni)
        S4 = self._S4(ksi, ni)
        S5 = self._S5(ksi, ni)
        S6 = self._S6(ksi, ni)

        S = np.hstack(([S1 * np.eye(2), S2 * np.eye(2), S3 * np.eye(2), S4 * np.eye(2), S5 * np.eye(2), S6 * np.eye(2)]))
        return S

    def _S1(self, ksi, ni):
        S1 = 1 - 3 * ksi**2 + 2 * ksi**3
        return S1

    def _S2(self, ksi, ni):
        S2 = self.L * (ksi - 2 * ksi**2 + ksi**3)
        return S2

    def _S3(self, ksi, ni):
        S3 = self.L * (ni - ksi * ni)
        return S3

    def _S4(self, ksi, ni):
        S4 = 3 * ksi**2 - 2 * ksi**3
        return S4

    def _S5(self, ksi, ni):
        S5 = self.L * (-ksi**2 + ksi**3)
        return S5

    def _S6(self, ksi, ni):
        S6 = self.L * ksi * ni
        return S6



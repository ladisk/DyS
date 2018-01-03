# coding=utf-8
"""

created by: lskrinjar
date of creation: 25/07/2016
time of creation: 21:40
"""
import time
from pprint import pprint
import itertools
import numpy as np
import vtk
import inspect


from MBD_system.body.body import Body
from MBD_system.q2R_i import q2R_i
from MBD_system.q2q_body import q2q_body
from MBD_system.gaussian_quadrature_on_arbitrary_intervals import change_interval


class Beam2D(Body):
    """
    classdocs
    """

    def __init__(self, id, node_id_i, node_id_j, mesh, nodes=None, parent=None):
        """

        """
        super(Beam2D, self).__init__(name="", parent=parent)

        #   pointer to mesh object
        self.mesh = mesh

        #   type of body
        self.body_type = "finite element"

        #   body id
        # self.body_id = self._count()

        #   element id
        self.element_id = id

        #   name
        self._name = "Beam2D_"+str(self.element_id)

        #   beam element nodes
        #   node ids
        self.node_id_i = self._node_id_i = node_id_i
        self.node_id_j = self._node_id_j = node_id_j

        #   node values
        if hasattr(self.mesh, "nodes0"):
            self.node_i = self.mesh.nodes0[node_id_i]
            self.node_j = self.mesh.nodes0[node_id_j]

        elif nodes is not None:
            self.node_i = self.R = nodes[node_id_i]
            self.node_j = nodes[node_id_j]

        else:
            raise ValueError, "Something is wrong!"

        #   list
        #   id list for mesh attribute "_nodes"
        self._node_id_list = [self._node_id_i, self._node_id_j]
        #   id list for mesh attribute "nodes" - this are ordered _nodes
        self.node_id_list = [self._node_id_i, self._node_id_j]

        #   parameter ksi list
        self._ksi_list = [0, 1]
        #   beam geometry
        #   length of the beam
        self.L_vector = self.node_j - self.node_i
        self.L = np.linalg.norm(self.L_vector)

        #   theta
        self.theta[2] = np.arctan2(self.L_vector[1], self.L_vector[0])

        #   color
        self.color = np.random.rand(3)

        #   geometry
        self.w = 0.
        self.h = 0.

        #   vector center of gravity
        self.r_g = None

        #   shape data (to display deformed shape)
        self.geometry_nodes = []
        #   number of subdivisions
        self.n_geometry_nodes = Beam2D.get_n_geometry_nodes()

        #   get geometry data from mesh
        self.module_of_elasticity = None
        self.density = None
        if self.mesh is not None:
            if hasattr(self.mesh._parent, "w"):
                self.w = self.mesh._parent.w

            if hasattr(self.mesh._parent, "h"):
                self.h = self.mesh._parent.h

            #   density
            if hasattr(self.mesh._parent, "density"):
                #   get density from flexible body object
                self.density = self.mesh._parent.density
            else:
                #   get density from mesh object
                self.density = self.mesh.density

            if hasattr(self.mesh._parent, "module_of_elasticity"):
                self.module_of_elasticity = self.mesh._parent.module_of_elasticity
            else:
                self.module_of_elasticity = self.mesh.module_of_elasticity

        #   evaluate volume of beam element
        if self.mesh is not None:
            if self.mesh.cross_section.area is None:
                self.mesh.cross_section.evaluate()

            #   volume
            self.volume = self.mesh.cross_section.area * self.L

            #   mass of beam element
            self.mass = self.density * self.volume

        #   predefined attributes for matrices
        self.B = None
        self.T = None
        self.J0 = None

        #   slope discontinuity
        self.g_1 = None
        #   transformation matrix
        self.A = None

        #   defined attributed, values are defined in subclasses
        self.e_n = None
        self.N = None
        #   integral boarders
        self.a = 0.
        self.b = 1.
        self.Izz = None

    @staticmethod
    def get_n_geometry_nodes():
        """

        :return:
        """
        return 41

    def _evaluate_S(self, ksi):
        """
        Evaluate shape function matrix (defined at every subclass)
        :return:
        """
        return np.empty([])

    def evaluate_r(self, x, e_i=None):
        """
        Evaluate shape function matrix (defined at every subclass)
        :return:
        """
        return np.zeros(2)

    def set_vtk_data(self):
        """


        :return:
        """
        #   vtk properties
        #   create a vtkPoints object and store the points in it
        self.nodes = vtk.vtkPoints()
        self.nodes.InsertNextPoint(np.append(self.mesh.nodes[self.node_id_i], 0.))
        self.nodes.InsertNextPoint(np.append(self.mesh.nodes[self.node_id_j], 0.))

        #   create the first line (between node i and j)
        self.element = vtk.vtkLine()
        self.element.GetPointIds().SetId(0, 0)  # the second 0 is the index of the Origin in the vtkPoints
        self.element.GetPointIds().SetId(1, 1)  # the second 1 is the index of P0 in the vtkPoints

        #   reate a cell array to store the lines in and add the lines to it
        self.line = vtk.vtkCellArray()
        self.line.InsertNextCell(self.element)

        #   create a polydata to store everything in
        self.vtk_linesPolyData = vtk.vtkPolyData()
        self.vtk_linesPolyData.SetPoints(self.nodes)
        self.vtk_linesPolyData.SetLines(self.line)

        self.vtk_mapper = vtk.vtkPolyDataMapper()
        self.vtk_mapper.SetInputData(self.vtk_linesPolyData)

        self.vtk_actor = vtk.vtkActor()
        self.vtk_actor.SetMapper(self.vtk_mapper)
        self.vtk_actor.GetProperty().SetLineWidth(10)
        self.vtk_actor.GetProperty().SetColor(self.color)

        if self.mesh._parent._visible_elements:
            self.vtk_actor.VisibilityOn()
        else:
            self.vtk_actor.VisibilityOff()

    def update_vtk_data(self, t, q):
        """

        :return:
        """
        #   absolute nodal coordinates of a beam element
        self.e = q2R_i(q, self.body_id)
        #   body coordinates R
        self.R[0:2] = self.e

        self.vtk_actor.SetPosition(self.R)

    def evaluate_geometry(self, e):
        """

        :param e:   vector of nodal coordinates of an element
        :return:    None
        """
        #   vector of nodes (x, y) on beam element
        x = np.zeros([self.n_geometry_nodes, 2])
        x[:, 0] = np.linspace(0, 1, num=self.n_geometry_nodes)

        #   predefine empty matrix
        self.geometry_nodes = np.zeros_like(x)

        for i in range(0, self.n_geometry_nodes):
            self.geometry_nodes[i, :] = np.dot(self._evaluate_S(x[i, 0]), e)

    def evaluate_e_i(self, e_b=None):
        """
        Function evaluates absolute nodal coordinates of a beam element in a mesh
        :param e_b:   vector of a mesh body
        :return:
        """
        # print "evaluate_e_i()"
        # curframe = inspect.currentframe()
        # calframe = inspect.getouterframes(curframe, 2)
        # print 'caller name:', calframe[1][3]
        #   vector of a mesh body
        if e_b is None:
            e_b = self.mesh.evaluate_q()
        # print "T ="
        # print self.evaluate_T().shape
        # print "B ="
        # print self.evaluate_B().shape
        # print "e_b =", e_b
        e_i = reduce(np.dot, [self.evaluate_T(), self.evaluate_B(), e_b])
        # print "e_i =", e_i
        return e_i

    def evaluate_CM(self, e_i):
        """

        :param e_i:
        :return:
        """
        m = np.array([[0.5, 0., self.L / 12., 0., 0.5, 0., -self.L / 12., 0.],
                      [0., 0.5, 0., self.L / 12., 0., 0.5, 0., -self.L / 12.]])

        CM = np.dot(m, e_i)

        return CM

    def evaluate_Q_g(self, g):
        """
        Function evaluates gravity force vector for finite element numerically with gauss quad
        :return:
        """
        #   gravity vector
        if len(g) == 3:
            g = g[0:2]

        #   predefine zero array
        Q_g = np.zeros(self.e_n)

        #   get beam absolute nodal coordinates
        # e_i = self.evaluate_e_i()

        #   get coordinates and weight of legendre polynomials for gauss quadrature
        x, w = np.polynomial.legendre.leggauss(self.N)

        for x_i, w_i in zip(x, w):
            x_i_interval, ab_i_interval = change_interval(x_i, self.a, self.b)
            Q_g_i = w_i * np.dot(g, self._evaluate_S(x_i_interval)) * ab_i_interval

            Q_g += Q_g_i

        return self.mass * Q_g

    def get_geometry_nodes(self, e=None):
        """

        :param e:   vector of absolute nodal coordinates of a mesh of a flexible body
        :return:
        """
        if e is None:
            if self.mesh is None:
                e = self.e0
            else:
                e = self.mesh.e



        # if len(e) != len(self.e0):
            #   vector of a flexible body
            # e_b = q2q_body(e, self.mesh._parent.body_id)

        #   vector of i-th beam element
        # print "self.evaluate_T() =", self.evaluate_T()
        # print "self.evaluate_B() =", self.evaluate_B()
        # print "e =", e
        if self.mesh is not None:
            e_i = reduce(np.dot, [self.evaluate_T(), self.evaluate_B(), e])

        else:
            e_i = e

        # else:
        #     if self.mesh is None:
        #         e_i = self.evaluate_q()
        #     else:
        #         e_i = self.mesh.evaluate_q()

        self.evaluate_geometry(e_i)

        return (self.geometry_nodes - self.geometry_nodes[0, :]) + self.geometry_nodes[0, :]

    def evaluate_M(self):
        """

        :return:
        """
        return None

    def evaluate_K(self, e):
        """

        :param e:
        :return:
        """

    def _evaluate_K_t(self):
        """
        Function evaluates stiffness matrix of a beam element - defined in subclass
        :return:
        """
        return None

    def _evaluate_K_l(self, e_i):
        """
        See ref: development of simple models for the elastic forces in the absolute nodal co-ordinate formulation (doi:10.1006/jsvi.1999.2935)
        Defined in subclass
        :param e:
        :return:
        """
        return None

    def evaluate_L(self, e_i):
        """
        Function evaluates length of the element with current vector of absolute nodal coordinates
        :return:
        """
        d = np.sqrt(((e_i[4] - e_i[0])**2) + (e_i[5] - e_i[1])**2)
        return d

    def _evaluate_dS_dksi(self, ksi):
        """
        This function is defined in subclass
        :return:
        """
        return None

    def _evaluate_ddS_ddksi(self, ksi):
        """
        This function is defined in subclass
        :return:
        """
        return None

    def evaluate_Q_s(self, q):
        """

        :param q:
        :return:
        """
        return None

    def evaluate_Q_e(self, q):
        """

        :param q:
        :return:
        """
        return None

    def evaluate_Q_e_M(self, q, M, ksi):
        """

        :param q:
        :param M:
        :param ksi:
        :return:
        """
        return None

    def evaluate_B(self):
        """
        Boolean matrix of the finite element always has a number of rows
        equal to the number of the finite element nodal coordinates and a number of
        columns equal to the number of the body nodal coordinates
        :return:
        """
        return None

    def evaluate_T(self):
        """
        Transformation matrix of a finite element relates the original element coordinates to the new set of coordinates of a body (mesh or structure).
        This matrix is not necessarily a square matrix.
        General Method for Modeling Slope Discontinuities and T-Sections Using ANCF Gradient Deficient Finite Elements
        :return:
        """
        return None

    def evaluate_C_q_fixed(self, node_id):
        """
        Overwritten in subclass
        :param node_id:
        :return:
        """
        return None

    def evaluate_C_q_hinged(self, node_id):
        """
        Overwritten in subclass
        :param node_id:
        :return:
        """
        return None

    def evaluate_kinetic_energy(self, q):
        """

        :return:
        """

    def evaluate_potential_energy(self, e_i, gravity):
        """

        :return:
        """
        self.evaluate_center_of_gravity(e_i)
        Ep = self.mass * np.dot(gravity[0:2], self.r_g)

        return Ep

    def evaluate_strain_energy(self, q):
        """

        :return:
        """

    def evaluate_center_of_gravity(self, e_i):
        """

        :param e_i: vector of absolute nodal coordinates of an element
        :return:
        """
        #   matrix
        Mg = np.array([[0.5, 0., self.L / 12., 0., 0.5, 0., -self.L / 12., 0.],
                       [0., 0.5, 0., self.L / 12., 0., 0.5, 0., -self.L / 12.]])

        self.r_g = np.dot(Mg, e_i)

        return self.r_g

    def plot(self, ax, color=None, show_id=True, linewidth=1):
        """

        :return:
        """
        if color is None:
            color = self.color

        xy_nodes = self.get_geometry_nodes()
        x = xy_nodes[:, 0]
        y = xy_nodes[:, 1]
        ax.plot(x, y, color=color, linewidth=linewidth)

        #   element id
        if show_id:
            ax.text(x[len(x)/2], y[len(y)/2], "ID: " + str(self.element_id), color=color)


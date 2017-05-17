"""

created by: lskrinjar
date of creation: 27/07/2016
time of creation: 16:28
"""
import itertools
import os
from pprint import pprint


import numpy as np
import scipy as sp
import vtk
from matplotlib import pyplot as plt


from MBD_system.MBD_system_items import BodyItem
import read_body_data_file.read_body_data_file as read_body_data_file
from MBD_system.Ai_ui_P import Ai_ui_P_vector
from MBD_system.extract_from_dictionary_by_string_in_key import extract_from_dictionary_by_string_in_key
from geometry.geometry import Geometry
from geometry.geometry_2D import Geometry2D
from MBD_system.q2dR_i import q2dR_i
from MBD_system.q2dtheta_i import q2dtheta_i
from MBD_system.q2R_i import q2R_i
from MBD_system.body.body import Body
from MBD_system.ancf.mesh.mesh import Mesh
from MBD_system.ancf.mesh.finite_element.beam_2d_euler_bernoulli import Beam2DEulerBernoulli
from MBD_system.ancf.mesh.finite_element.beam_2d_shear_deformable import Beam2DShearDeformable
from global_variables import GlobalVariables
from MBD_system.q2q_body import q2q_body
from MBD_system.q2dq_body import q2dq_body
from MBD_system import transform_cs_type


class FlexibleBody(Body):
    """
    classdocs
    """

    def __init__(self, name, mesh_file_path=None, file_path=None, _dict={}, cross_section_type=None, parent=None):
        """
        Constructor of body class
        :param name:                    body name (string)
        :param filename:                absolute path file of body properties
        :param R:                       a vector of positions R (numpy array) in m
        :param theta:                   orientation angles (numpy array) in degrees
        :param dR:                      a vector of velocities (numpy array) in m/s
        :param dtheta:                  a vector of angular velocities (numpy array) in deg/s
        :param color:                   a color vector (RGB)
        :param transparent:             transparency (float) a value from 0 to 1
        :param visible:                 true or false
        :param display_style:           a display style (as string) options: filled, wireframe, points-check last option
        :param properties_file:         a path to mass and geometry properties data in .dat file (todo)
        :param geometry_data_file:      a path to geometry .stl or .obj file (todo)
        """
        super(FlexibleBody, self).__init__(name=name, file_path=file_path, parent=parent)

        #   body id
        self.body_id = self._count()

        #   body type
        self.body_type = "flexible body"
        self.cross_section_type = cross_section_type

        #   visualization properties
        self.color = np.random.rand(3)
        self._visible_elements = False #False or True
        # self._element_visibility()

        #   predefine mesh attributes
        self.mesh_filename = mesh_file_path
        #   folder name to store deformed shaped at every time step (just for testing)
        self._results_folder = self._name + "_results"

        #   flexible body additional attributes
        self.K = None
        self.dissipation_coefficient = None

        #   energy
        self._elastic_strain_energy = 0.

        #   mesh object as attribute
        self.mesh = Mesh(parent=self)

        #    read body properties file
        if self.file_path is not None:
            if os.path.isfile(self.file_path):
                if self._dict == {}:
                    self._dict = read_body_data_file.read_body_data_file(self.file_path)

                self.add_attributes_from_dict(self._dict)

            else:
                raise ValueError, "File not found!"

        else:
            print "File path not correct! File path is: %s"%self.file_path

        #   properties of cross section of a mesh
        self.mesh.cross_section.type = self.cross_section_type
        if self.mesh.cross_section.type is not None:
            self.mesh.cross_section.evaluate()

        #   body coordinates
        self.q_i_size = None
        self.q_i_dim = None

        #   check if mesh file exists and read it
        if self.mesh_filename is not None:
            #   absolute file path
            if self._parent is not None:
                file_path = os.path.join(self._parent._parent.MBD_folder_abs_path, self.mesh_filename)
                self._results_folder_path = os.path.join(self._parent._parent.MBD_folder_abs_path, self._results_folder)
            else:
                file_path = os.path.join(os.getcwd(), self.mesh_filename)
                self._results_folder_path = os.path.join(os.getcwd(), self._results_folder)

            if os.path.isfile(file_path):
                self.mesh.read_file(file_abs_path=file_path)
            else:
                print Warning, "Filename %s does not exist! Check filename!"%file_path

            #   size of vector q for a flexible body - mesh
            self.q_i_dim = self.mesh.n_NC
            self.e_n = self.mesh.n_NC
        else:
            print Warning, "Attribute mesh_filename is None! Mesh file not read!"

        #   transform data to ANC in GCS
        self.transform2ANC()

        #   check to save results option
        self._save_results(self._results_folder_path)

    def _excel_header_q(self):
        """

        :return:
        """
        id = str(self.body_id)

        hdr = []
        for i, (node, node_dim) in enumerate(zip(self.mesh.nodes, self.mesh.node_dim)):
            if node_dim == 4:
                hdr_i = [id + "_e" + str(i) + "_x",
                         id + "_e" + str(i) + "_y",
                         id + "_ex" + str(i) + "_x",
                         id + "_ex" + str(i) + "_y"]

            if node_dim == 6:
                hdr_i = [id + "_e" + str(i) + "_x",
                         id + "_e" + str(i) + "_y",
                         id + "_eX" + str(i) + "_x",
                         id + "_eX" + str(i) + "_y",
                         id + "_eY" + str(i) + "_x",
                         id + "_eY" + str(i) + "_y"
                         ]

            hdr += hdr_i

        return hdr

    def _excel_header_dq(self):
        """

        :return:
        """
        id = str(self.body_id)

        hdr = []
        for i, (node, node_dim) in enumerate(zip(self.mesh.nodes, self.mesh.node_dim)):
            if node_dim == 4:
                hdr_i = [id + "_de" + str(i) + "_x",
                         id + "_de" + str(i) + "_y",
                         id + "_dex" + str(i) + "_x",
                         id + "_dex" + str(i) + "_y"]

            if node_dim == 6:
                hdr_i = [id + "_de" + str(i) + "_x",
                         id + "_de" + str(i) + "_y",
                         id + "_deX" + str(i) + "_x",
                         id + "_deX" + str(i) + "_y",
                         id + "_deY" + str(i) + "_x",
                         id + "_deY" + str(i) + "_y"
                         ]

            hdr += hdr_i

        return hdr

    def _element_visibility(self):
        """
        Method changes visibility of finite elements in a mesh. Each undeformed element has different color.
        :return:
        """
        if self._visible_elements:
            self._visible_elements = False
        else:
            self._visible_elements = True

        for element in self.mesh.elements:
            if element.vtk_actor.GetVisibility():
                element.vtk_actor.VisibilityOff()
            else:
                element.vtk_actor.VisibilityOn()

    def set_vtk_data(self):
        """

        :return:
        """
        #   evaluate mesh geometry at initial undeformed state
        self.mesh.evaluate_geometry_nodes()

        #   points
        self.geometry_nodes_vtk = vtk.vtkPoints()
        self.geometry_nodes_vtk.SetNumberOfPoints(len(self.mesh.geometry_nodes))

        #   vertices
        self.vtkCells = vtk.vtkCellArray()

        #   add nodes
        for i, node in enumerate(self.mesh.geometry_nodes):
            self.geometry_nodes_vtk.SetPoint(i, node[0], node[1], 0.0)

        #   vtkCellArray is a supporting object that explicitly represents cell connectivity.
        #   The cell array structure is a raw integer list of the form:
        #   (n,id1,id2,...,idn, n,id1,id2,...,idn, ...) where n is the number of points in
        #   the cell, and id is a zero-offset index into an associated point list.
        lines = vtk.vtkCellArray()
        lines.InsertNextCell(len(self.mesh.geometry_nodes))
        for i in xrange(0, len(self.mesh.geometry_nodes)):
            lines.InsertCellPoint(i)

        #   vtkPolyData is a data object that is a concrete implementation of vtkDataSet.
        #   vtkPolyData represents a geometric structure consisting of vertices, lines, polygons, and/or triangle strips
        self.vtkPolyData = vtk.vtkPolyData()
        self.vtkPolyData.SetPoints(self.geometry_nodes_vtk)
        self.vtkPolyData.SetLines(lines)
        self.vtkPolyData.SetVerts(self.vtkCells)

        self.polygonMapper = vtk.vtkPolyDataMapper()
        self.polygonMapper.SetInputData(self.vtkPolyData)
        self.polygonMapper.Update()

        #   vtk actor
        self.vtk_actor = vtk.vtkActor()
        self.vtk_actor.SetMapper(self.polygonMapper)
        #   set color
        self.vtk_actor.GetProperty().SetColor(self.color)

        # self.vtk_actor.SetPosition(self.R)
        # self.vtk_actor.SetOrientation(np.rad2deg(self.theta))

        #   marker
        self.set_vtk_LCS()

    def _update_vtk_data(self, t, e):
        """

        :return:
        """
        #   absolute nodal coordinates of a flexible body (mesh, structure)
        self.e_b = self.q = e[np.sum(GlobalVariables.q_i_dim[0:self.body_id]):np.sum(GlobalVariables.q_i_dim[0:self.body_id + 1])]
        # self.e_b = self.q =
        #   evaluate mesh geometry at current state
        self.mesh.evaluate_geometry_nodes(e_b=self.e_b)

        # if self.save_results:
        #     fig = plt.figure(figsize=(6, 5),
        #                      dpi=100,
        #                      facecolor='w',
        #                      edgecolor='k')
        #
        #     dxy = 1.1 * abs(np.array(self.mesh.geometry_nodes)).max()
        #     plt.xlim([-dxy, dxy])
        #     plt.ylim([-dxy, dxy])
        #     plt.plot(np.array(self.mesh.geometry_nodes)[:, 0], np.array(self.mesh.geometry_nodes)[:, 1])
        #
        #     #   simulation step
        #     plt.figtext(0.1, 0.1, "step="+str(self.step).zfill(4), fontdict=None)
        #     plt.grid(True)
        #     filename = "step_" + str(self.step).zfill(4) + ".png"
        #     # filename = "step_" + str(self.step).zfill(4) + ".png"
        #     fig.savefig(os.path.join(self._results_folder_path, filename))
        #     fig.clf()

        #   replace nodes
        for i, node in enumerate(self.mesh.geometry_nodes):
            node = np.append(node, 0.)
            #   set point
            self.geometry_nodes_vtk.SetPoint(i, node)

        #   data is changed
        self.geometry_nodes_vtk.Modified()
        self.polygonMapper.Update()

        #   marker
        # R = np.append(self.e_b[0:2], 0.) + self.R
        # theta = np.array([0., 0., transform_cs_type.transform_cartesian2polar(self.e_b[2], self.e_b[3])[1]])
        # self.update_vtk_LCS(R=R, theta=theta)
        self.R = self.mesh.geometry_nodes[0,:]
        # self.update_vtk_LCS(R=self.mesh.geometry_nodes[0,:])

        #   update nodes
        # self.mesh.update_vtk_data(self.e_b)

    def update_coordinates_and_angles_2D(self, q_i):
        """
        Only coordinates are updated as no angles are used in ANCF formulation for flexible body
        :return:
        """
        self.R[0:2] = q_i[0:2]

    def transform2ANC(self):
        """
        Function transforms nodes from mesh coordinates to ANC in GCS
        :return:
        """
        if (self.R != np.zeros(3)).any() or (self.theta != np.zeros(3)).any():
            self.mesh.transform2ANC(self.R[0:2], self.theta[2])

    def evaluate_q_i_size(self):
        """

        :return:
        """
        self.q_i_size = len(self.mesh.evaluate_q())
        return self.q_i_size

    def evaluate_q0(self):
        """

        :return:
        """
        self.q = self.mesh.evaluate_q()
        return self.q

    def evaluate_q(self):
        """

        :rtype: object
        :return:
        """
        self.q = self.mesh.evaluate_q()
        return self.q

    def evaluate_dq(self):
        """

        :return:
        """
        self.dq = np.zeros_like(self.q)
        return self.dq

    def evaluate_dq0(self):
        """

        :return:
        """
        self.dq = self.mesh.evaluate_dq0()
        return self.dq

    def evaluate_M_size(self):
        """

        :return:
        """
        self.M_size = len(self.mesh.evaluate_q())

        return self.M_size

    def evaluate_M(self):
        """

        :return:
        """
        # print "mass =", self.mesh.evaluate_mass()
        # print "number of elements =", self.mesh.n_e
        # print "number of element nodal coordinates =", self.mesh.elements[0].e_n
        # print "number of nodes =", len(self.mesh.nodes0)
        # print "number of nodal coordinates in a mesh =", self.mesh.n_NC
        self.e_n = self.mesh.n_NC
        self.M = self.mesh.evaluate_M()
        # print "size =", np.shape(self.M)
        return self.M

    def evaluate_CM(self, q):
        """

        :param q:
        :return:
        """
        e_i = q2q_body(q, self.body_id)

        self.CM = self.mesh.evaluate_CM(e_i)
        return self.CM

    def evaluate_Q_e(self, q):
        """
        Function evaluates vector of generalized external forces on flexible body
        :param q:
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

    def evaluate_r(self, q, node_id=None, element_id=None, ksi=None):
        """

        :param q:       vector of absolute nodal coordinates
        :param node_id: id of node
        :return:
        """
        #   vector of body nodal coordinates
        e_b = self.evaluate_q2q_body(q)

        #   vector of coordinates of a node in GCS
        #   body i
        #   node j
        rP_ij = self.mesh.evaluate_r(e_b, node_id, element_id=element_id, ksi_element=ksi)
        return rP_ij

    def evaluate_dr(self, q, element_id=None, ksi=None):
        """

        :param q:
        :param node_id:
        :param element_id:
        :param ksi:
        :return:
        """
        #   vector of body nodal coordinates
        de_b = q2dq_body(q, self.body_id)
        #   vector of velocity at position ksi on deformable body i in GCS
        drP_i = self.mesh.evaluate_dr(de_b, element_id=element_id, ksi_element=ksi)
        return drP_i

    def e_node(self, q, node_id):
        """

        :param node_id:
        :return:
        """
        e_i = q2q_body(q, self.body_id)

        e_node = self.mesh.e_node(e_i, node_id)
        return e_node

    def evaluate_C_q_fixed(self, node_id):
        """

        :param node_id:
        :return:
        """
        #   constrain on absolute nodal coordinates e1, e2 (x, y)
        #   find element that has selected node_id inside
        element = [element for element in self.mesh.elements if node_id in element.node_id_list][0]

        return element.evaluate_C_q_fixed(node_id)

    def evaluate_C_q_hinged(self, node_id):
        """

        :param node_id:
        :return:
        """
        #   constrain on absolute nodal coordinates e1, e2 (x, y)
        #   find element that has selected node_id inside
        # element = [element for element in self.mesh.elements if node_id in element.node_id_list][0]

        #   predefine empty matrix C_q
        C_q = np.zeros([2, self.q_i_dim])

        #   ones at nodal coordinates of position
        C_q[:, sum(self.mesh.node_dim[0:node_id]):sum(self.mesh.node_dim[0:node_id]) + 2] = np.eye(2)

        return C_q

    def evaluate_Q_g(self, q=None, gravity=None):
        """
        Function evaluates vector of generalized distributed gravity forces of a flexible body
        :return:
        """

        if gravity is None:
            gravity = self._parent._parent.gravity
        else:
            gravity = gravity[0:2]

        if q is None:
            q_i = self.evaluate_q()
        else:
            q_i = q2q_body(q, self.body_id)

        self.Q_g = self.mesh.evaluate_Q_g(g=gravity, q=q_i)

        return self.Q_g

    def print_Q_g(self):
        """

        :return:
        """
        Q_g = self.evaluate_Q_g()
        print "Q_g ="
        print Q_g

    def evaluate_K(self, e=None):
        """

        :return:
        """
        #   vector of ANC of a flexible body
        if e is not None:
            e = q2q_body(e, self.body_id)
        else:
            e = self.mesh.evaluate_q()

        self.K = self.mesh.evaluate_K(e)
        return self.K

    def print_K(self):
        """

        :return:
        """
        print "Stiffness matrix of %s"%self._name
        print self.evaluate_K()

    def evaluate_Q_s(self, q):
        """
        Function evaluates elastic strain forces of flexible body (mesh, structure)
        :param q_i:
        :return:
        """
        #   vector of absolute nodal coordinates of a flexible body
        e_i = q2q_body(q, self.body_id)
        # Q_s_i = np.zeros_like(e_i)

        #   stiffness matrix
        # self.K = self.evaluate_K(e_i)

        #   elastic strain forces of body flexible body i
        # Q_s_i = np.dot(self.K, e_i)

        Q_s_i = self.mesh.evaluate_Q_s(e_i)

        return Q_s_i

    def evaluate_Q_dis(self, q):
        """

        :param q:
        :return:
        """
        #   vector of absolute nodal coordinates of a flexible body
        de_i = q2dq_body(q, self.body_id)

        #   dissipation of body flexible body i
        Q_dis = self.dissipation_coefficient * np.dot(self.M, de_i)

        return Q_dis

    def check_slope_discontinuity(self, parent):
        """

        :param parent:
        :return:
        """
        self.mesh.check_slope_discontinuity(parent)

    def print_nodes(self):
        """

        :return:
        """
        print "test1"
        e = self.evaluate_q()
        print "Geometry nodes:"
        print self.mesh.evaluate_geometry_nodes(e=e)

    def evaluate_mesh_data_at_position(self, L):
        """

        :return:
        """
        ksi, element_id, node_id = self.mesh.evaluate_ksi_element_id_node_id_at_position(L)

        return ksi, element_id, node_id

    def evaluate_kinetic_energy(self, q):
        """

        :return:
        """
        dq_b = self.evaluate_q2dq_body(q)

        self._kinetic_energy = self.mesh.evaluate_kinetic_energy(dq_b)

        return self._kinetic_energy

    def evaluate_potential_energy(self, q, gravity):
        """

        :return:
        """
        q_b = self.evaluate_q2q_body(q)

        self._potential_energy = self.mesh.evaluate_potential_energy(q_b, gravity)

        return self._potential_energy

    def evaluate_elastic_strain_energy(self, q, q_i_dim=None):
        """

        :return:
        """
        # print "q =", type(q)
        # e = q2q_body(q, self.body_id, q_i_dim=q_i_dim)

        e = self.evaluate_q2q_body(q)

        self._elastic_strain_energy = self.mesh.evaluate_elastic_strain_energy(e)

        return self._elastic_strain_energy


if __name__ == "__main__":
    # filename = "mesh_2d_2elements.mesh"
    # filename = "frame_mesh_2d_2elements.mesh"
    filename = "moving_contact_line_2elements.mesh"
    body = FlexibleBody("test_mesh")
    body.mesh.element_type = "beam 2D euler-bernoulli"
    #   assign geometry and mass properties
    body.mesh.cross_section.type = "rectangle"
    body.mesh.cross_section.B = .1
    body.mesh.cross_section.H = .5

    body.mesh.density = 7800.
    body.mesh.module_of_elasticity = 2.07E+11

    body.mesh.read_file(file_abs_path=filename)

    L_position = 18.5E-3
    print body.evaluate_mesh_data_at_position(L_position)

    # print "e ="
    # print body.evaluate_q()
    #
    # for element in body.mesh.elements:
    #     print element.evaluate_e_i()





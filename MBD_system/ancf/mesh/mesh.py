"""

created by: lskrinjar
date of creation: 25/07/2016
time of creation: 22:34
"""
import os
from pprint import pprint
import numpy as np
import scipy as sp
from matplotlib import pyplot as plt
import matplotlib as mpl
import inspect
import vtk


from finite_element.beam_2d_euler_bernoulli import Beam2DEulerBernoulli
from finite_element.beam_2d_shear_deformable import Beam2DShearDeformable
from mesh_options import MeshOptions
from cross_section.cross_section import CrossSection
from MBD_system.q2q_body import q2q_body
from MBD_system.ancf.mesh.slope_discontinuity import SlopeDiscontinuity
from MBD_system.transform_cs import uP_lcs2gcs


class Mesh(object):
    """
    Mesh object constructor
    """

    def __init__(self, element_type=None, parent=None):
        """

        :return:
        """
        #   parent
        self._parent = parent

        #   geometry
        self.lines = []

        #   mesh options
        self.mesh_options = MeshOptions(parent=self)

        #   cross section
        self.cross_section = CrossSection(parent=self)

        #   material properties (defined object attributes wih predefined values)
        self.density = None
        self.module_of_elasticity = None
        self.poisson_ratio = None

        #   matrices
        self.M = None
        self.K = None
        self.K_constant_matrix = False

        #   energy
        self._kinetic_energy = 0.
        self._potential_energy = 0.
        self._strain_energy = 0.

        #   time integration (dependent) properties
        self.step = None

        #   scale data
        self.scale = 1E-3

        #   initial geometry data (unsorted) from gmsh file
        self.nodes0 = []
        #   list of node ids
        self.node_id_list0 = []
        self.normals0 = []

        #   list of mesh nodes (ordered)
        self.nodes = []
        #   list of nodes from which a mesh is constructed, duplicate nodes are removed automatically during construction
        self.node_id_list = []
        self.normals = []

        #   total length of beam elements
        self.L0 = 0.

        #   list of finite elements
        self.elements = []

        #   list of nodal coordinates at every node in a mesh
        self.node_dim = []

        #   type of elements
        self.element_type = element_type
        if hasattr(self._parent, "element_type"):
            self.element_type = self._parent.element_type

            if self.element_type == Beam2DEulerBernoulli.get_element_type():
                self.element_e_n = Beam2DEulerBernoulli.get_e_size()

            elif self.element_type == Beam2DShearDeformable.get_element_type():
                self.element_e_n = Beam2DShearDeformable.get_e_size()

            else:
                print Warning, "Element type not defined!"

            self.node_e_n = self.element_e_n / 2

        else:
            self.element_e_n = None

        #   number of nodes in a mesh
        self.n_n = None

        #   number of finite elements in a mesh
        self.n_e = len(self.elements)

        #   vector of nodal coordinates of all elements before element connection
        self.e_b = []

        #   vector of nodal coordinates of a mesh
        self.e = []

        #   number of nodal coordinates
        self.n_NC = 0

        #   predefined empty list of slope discontinuities
        self.slope_discontinuities = []

        self.n_geometry_nodes = 0

        #   filename of a mesh file
        self.file_abs_path = None

        #   visualization properties
        self.vtk_actor = None
        self.vtk_mapper = None
        self.vtk_poly_data = None
        self.vtk_cells = None
        self.vtk_points = None

        self._node_size = 2

    def scale(self, scale=None):
        """
        Function scales geometry
        :param scale:
        :return:
        """
        if scale is not None:
            self.scale = scale

        self.nodes = self.nodes * self.scale

    def set_vtk_data(self):
        """
        Display mesh nodes
        :return:
        """
        self.vtk_poly_data = vtk.vtkPolyData()

        self.vtk_points = vtk.vtkPoints()

        self.vtk_mapper = vtk.vtkPolyDataMapper()
        self.vtk_mapper.SetInputData(self.vtk_poly_data)

        self.vtk_actor = vtk.vtkActor()
        self.vtk_actor.SetMapper(self.vtk_mapper)
        self.vtk_actor.GetProperty().SetLineWidth(10)

        self.vtk_cells = vtk.vtkCellArray()

        self.vtk_poly_data.SetPoints(self.vtk_points)
        # self.vtk_poly_data.SetVerts(self.vtk_cells)

        # for i, node in enumerate(self.nodes):
        #     node = np.append(node, 0.)
        #     pointId = self.vtk_points.InsertNextPoint(node)
        #     self.vtk_cells.InsertNextCell(i)
        #     self.vtk_cells.InsertCellPoint(pointId)

        self.vtk_actor.GetProperty().SetPointSize(self._node_size)

    def update_vtk_data(self, e):
        """

        :return:
        """
        #   update nodes
        # for i, node in enumerate(self.nodes):
        #     node = e[sum(self.node_dim[0:i]):sum(self.node_dim[0:i]) + 2]
        #     node = np.append(node, 0.)
        #     #   set point
        #     self.vtk_points.SetPoint(i, node)
        #
        # #   data is changed
        # self.vtk_points.Modified()
        # self.vtk_mapper.Update()

    def create_mesh(self):
        """
        Function creates mesh - list of elements (without reading mesh file)
        :return:
        """
        self.nodes, self.elements = self._create_line_mesh()

        #   number of finite elements in a mesh
        self.n_e = len(self.elements)

    def add_lines(self, lines):
        """

        :param lines:
        :return:
        """
        self.lines = lines

    def _create_line_mesh(self):
        """

        :return:
        """
        #   construct list of elements
        elements = []
        #   construct nodes
        self._nodes0 = []

        #   for every line in a list
        for line_indx, line in enumerate(self.lines):
            #   divide line by number of elements
            if self.mesh_options.number_of_elements is not None and self.mesh_options.element_size is None:
                pass

            #   divide line by element size
            elif self.mesh_options.number_of_elements is None and self.mesh_options.element_size is not None:
                #   assign property to mesh options
                self.mesh_options.number_of_elements = int(np.rint(line.length / self.mesh_options.element_size))

            else:
                raise ValueError, "Division properties are not defined!"

            #   a vector of divisions on a line
            dL = np.linspace(0., 1., num=self.mesh_options.number_of_elements + 1) * line.length

            for i in range(0, len(dL)):
                node = line.node_i + (line.tangent_unit * dL[i])
                if not self._nodes0:
                    self._nodes0.append(node)
                else:
                    #   check if node is already in a list and if not, append it
                    if np.all(np.any(np.array(self._nodes0) - node, axis=1)):
                        self._nodes0.append(node)

            for i in range(len(self._nodes0) - self.mesh_options.number_of_elements - 1, len(self._nodes0) - 1):
                node_id_i = i
                node_id_j = i + 1

                #   create finite element object
                if self.element_type == Beam2DEulerBernoulli.get_element_type():
                    id = len(elements)
                    element = Beam2DEulerBernoulli(id, node_id_i, node_id_j, self)

                elif self.element_type == Beam2DShearDeformable.get_element_type():
                    element = Beam2DShearDeformable()

                else:
                    raise ValueError, "Element type not corrent!"

                elements.append(element)

        return self._nodes0, elements

    def read_file(self, file_abs_path=None):
        """
        Read .mesh file
        :param      file_abs_path: absolute path to .mesh file
        :return:    None
        """
        if file_abs_path is not None:
            self.file_abs_path = os.path.abspath(file_abs_path)

        #   load data from file if already exists
        _path, self._file = os.path.split(self.file_abs_path)
        self._name, self._filetype = os.path.splitext(self._file)

        if self._filetype == ".mesh":
            self._nodes0, self.nodes, self.node_id_list, self.elements = self._read_mesh_file()

        #   number of nodal coordinates
        if self.element_type == Beam2DEulerBernoulli.get_element_type():
            self.n_NC = len(self.nodes0) * 4

        elif self.element_type == Beam2DShearDeformable.get_element_type():
            self.n_NC = len(self.nodes0) * 6

        else:
            print Warning, "Element type not supported!"

        #   number of finite elements in a mesh
        if self.elements:
            self.n_e = len(self.elements)
            #   number of subdivision nodes for deformed mesh
            self.n_geometry_nodes = self.n_e * self.elements[0].n_geometry_nodes
        else:
            print Warning, "Elements list is empty!"

        #   vector of all nodal coordinates before element connection
        for element in self.elements:
            self.e_b.append(element.e0)

        if self.elements:
            self.e_b = np.concatenate(self.e_b, axis=0)

        #   vector of nodal coordinates of a mesh
        self.e = self.evaluate_q()

        #   mesh geometry
        self.nodes, self._nodes0 = self._construct_nodes_2D(self.nodes)

    def _read_mesh_file(self, filename=None):
        """

        :param filename:
        :return:
        """
        if filename is None:
            filename = self.file_abs_path

        #   predefine list to save data in file
        with open(filename,'r') as _file:
            out = []
            for line in _file:
                # print line
                out.append(line.strip())
        _file.close()

        #   find vertices and create them
        if 'Vertices' in out:
            self.nodes0, self.node_id_list = self._read_mesh_file_nodes(out)

        else:
            raise ValueError('ReadMesh: no vertices defined.')

        #   find elements and create them
        if 'Edges' in out:
            elements, nodes, node_id_list = self._read_mesh_file_elements(out)
        else:
            raise ValueError('ReadMesh: no elements defined.')

        #   del variable when finished creating mesh properties
        del out

        print "Mesh read successfully from file: ", filename

        return self.nodes0, nodes, node_id_list, elements

    def _read_mesh_file_nodes(self, data):
        """

        :param data:
        :return:
        """
        ind = data.index('Vertices')+1
        self.n_n = int(data[ind].strip())

        nodes = []
        indx = []
        for i, node in enumerate(data[ind+1:ind+self.n_n + 1]):
            indx.append(i)
            _node = np.array(np.fromstring(node, dtype=float, count=-1, sep=' ')[0:3]) * self.scale
            nodes.append(_node)

        return nodes, indx

    def _read_mesh_file_elements(self, data):
        """
        Method reads mesh file and creates list of elements
        :param data:
        :return:
        """
        #   nodes list ordered as used in each element
        nodes = []
        node_id_list = []

        #   predefined list to store all elements in a mesh
        elements = []
        ind = data.index('Edges') + 1
        self.n_e = int(data[ind][0])

        #   read mesh file
        #   list of elements
        for i, element_str in enumerate(data[ind+1:ind+self.n_n]):
            #   read each row
            [node_id_i, node_id_j, unknown_param] = np.fromstring(element_str, dtype=int, count=-1, sep=' ')

            if self.element_type == Beam2DEulerBernoulli.get_element_type():
                #   node ids are shifted to start with 0 and not with 1 as in gmsh file
                node_id_i = node_id_i - 1
                node_id_j = node_id_j - 1

                #   create element object and append it to list
                element = Beam2DEulerBernoulli(i, node_id_i, node_id_j, self)
                elements.append(element)

        #   sort nodes by element id
        nodes, node_id_list = self.sort_nodes(elements)

        return elements, nodes, node_id_list

    def sort_nodes(self, elements):
        """
        Sort nodes by element.
        :return:
        """
        nodes = []
        node_id_list = []

        for i_element, element in enumerate(elements):
            for i_node, (node_id, node_id_str) in enumerate(zip(element.node_id_list, ["node_id_i", "node_id_j"])):
                node = self.nodes0[node_id]
                #   if list is empty
                if not nodes:
                    _id = len(nodes)
                    nodes.append(node)
                    node_id_list.append(_id)

                    setattr(element, node_id_str, _id)
                    element.node_id_list[0] = _id

                #   if list is not empty (after first node is appended to list)
                else:
                    if not np.equal(node, nodes).all(axis=1).any():
                        _id = len(nodes)
                        nodes.append(node)
                        node_id_list.append(_id)

                        setattr(element, node_id_str, _id)
                        element.node_id_list[i_node] = _id
                    else:
                        _id = np.where(np.all(nodes==node, axis=1))[0][0]

                        setattr(element, node_id_str, _id)
                        element.node_id_list[i_node] = _id

        return nodes, node_id_list

    def _construct_nodes_2D(self, nodes):
        """
        Function creates nodes attributes of mesh object
        :argument: nodes list of node items, node item is type numpy array
        :return:
        """
        nodes_2D = []
        for node in nodes:
            nodes_2D.append(node[0:2])

        return nodes_2D, nodes_2D

    def transform2ANC(self, R, theta):
        """
        Transform nodes
        :return:
        """
        for i, (node, node_id) in enumerate(zip(self.nodes, self.node_id_list)):
            self._nodes0[node_id] = self.nodes0[node_id] = self.nodes[node_id] = R + uP_lcs2gcs(u_P=self.nodes[node_id], theta=theta)

        for element in self.elements:
            element.theta[2] = theta

        self.evaluate_q()

    def evaluate_M(self):
        """

        :return:
        """
        #   if attribute is an empty list, run evaluate_q()
        if not self.node_dim:
            self.evaluate_q()

        M = np.zeros([self.n_NC, self.n_NC])
        for element in self.elements:
            M += reduce(np.dot, [element.evaluate_B().T, element.evaluate_T().T, element.evaluate_M(), element.T, element.B])

        self.M = M
        return self.M

    def evaluate_CM(self, e_b):
        """

        :param q:
        :return:
        """
        CM = np.zeros(2, dtype=float)
        M = 0.
        for element in self.elements:
            e_i = element.evaluate_e_i(e_b=e_b)

            CM_i = element.evaluate_CM(e_i) * element.mass
            CM += CM_i

            M += element.mass

        CM = CM / M

        return CM

    def evaluate_Q_s(self, e):
        """

        :param e:
        :return:
        """
        # print "evaluate_Q_s() @", __name__
        if self.K_constant_matrix:
            Q_s = np.dot(self.K, e)
        else:
            Q_s = np.zeros_like(e)

            for i, element in enumerate(self.elements):
                # print i
                e_i = element.evaluate_e_i(e_b=e)
                Q_s_i = element.evaluate_Q_s(e_i)
                Q_s += reduce(np.dot, [element.evaluate_B().T, element.evaluate_T().T, Q_s_i])

        # print "Q_s =", Q_s
        return Q_s

    def evaluate_mass(self):
        """

        :return:
        """
        mass = 0.

        for element in self.elements:
            mass += element.mass

        return mass

    def evaluate_K(self, e):
        """
        :param e:   vector of absolute nodal coordinates of a mesh (body)
        :return:
        """
        # print "evaluate_K()"
        # print "e =", len(e)

        K = np.zeros([self.n_NC, self.n_NC])

        for element in self.elements:
            e_i = element.evaluate_e_i(e_b=e)
            K_e = element.evaluate_K(e_i)
            K_e_augmented = reduce(np.dot, [element.evaluate_B().T, element.evaluate_T().T, K_e, element.T, element.B])
            K += K_e_augmented

        self.K = K
        return self.K

    def _construct_constant_matrix_K_t(self):
        """

        :return:
        """
        K = np.zeros([self.n_NC, self.n_NC])

        for element in self.elements:
            K_e = element._construct_constant_matrix_K_t()
            K_e_augmented = reduce(np.dot, [element.evaluate_B().T, element.evaluate_T().T, K_e, element.T, element.B])
            K += K_e_augmented

        self.K = K
        return self.K

    def _evaluate_K_t(self):
        """

        :return:
        """
        return None

    def _evaluate_K_l(self, e):
        """
        See ref: development of simple models for the elastic forces in the absolute nodal co-ordinate formulation (doi:10.1006/jsvi.1999.2935)
        :param e:
        :return:
        """
        return None

    def evaluate_Q_g(self, g, q=None):
        """
        :param q:   vector of body (mesh, structure) absolute nodal coordinates
        :param g:   gravity vector (2D (x, y))
        :return:
        """
        if q is None:
            q = self.evaluate_q()

        #   predefine zero array
        Q_g = np.zeros(self.n_NC)

        for element in self.elements:
            # Q_g_i = reduce(np.dot, [element.evaluate_T(), element.evaluate_B(), element.evaluate_Q_g(g)])
            Q_g_i = element.evaluate_Q_g(g)
            # print "Q_g_i =", Q_g_i
            #   augmented to mesh coonrdinates
            Q_g_i_augmented = reduce(np.dot, [element.evaluate_B().T, element.evaluate_T().T, Q_g_i])
            Q_g += Q_g_i_augmented
            # Q_g += np.dot(element.evaluate_B().T, element.evaluate_Q_g(g))

        # print "Q_g =", len(Q_g), Q_g
        return Q_g

    def evaluate_geometry_nodes(self, e_b=None):
        """
        Method creates geometry nodes of undeformed and deformed geometry during simulation
        :param e:       vector of body (mesh, structure) absolute nodal coordinates
        :param e_b:     gravity vector (2D (x, y))
        :return:
        """
        #   nodes
        #   x component: self.nodes[0,:] - col 0
        #   y component: self.nodes[1,:] - col 1
        # print "self.n_geometry_nodes =", self.n_geometry_nodes
        self.geometry_nodes = np.zeros([self.n_geometry_nodes, 2])
        for i, element in enumerate(self.elements):
            n_e = element.get_geometry_nodes(e=e_b)
            # print "n_e =", n_e
            self.geometry_nodes[i*element.n_geometry_nodes:(i + 1)*element.n_geometry_nodes, :] = n_e

        return self.geometry_nodes

    def evaluate_q(self):
        """
        Function constructs vector e(q) of finite element mesh and vector of absolute nodal coordinates at every node.
        If there is slope discontinuity at node additional slope vectors in structure (mesh) coordinates are defined.
        :return:
        """
        # print "-------------------------------------"
        # print "evaluate_q()", __name__
        curframe = inspect.currentframe()
        calframe = inspect.getouterframes(curframe, 2)
        # print 'caller name:', calframe[1][3]
        # print "self.elements =", self.elements
        # print "self.nodes =", self.nodes
        # print "self._nodes0 =", self._nodes0

        if self.element_e_n is None:
            if self.elements:
                if self.elements[0].finite_element_type == Beam2DEulerBernoulli.get_element_type():
                    self.element_e_n = Beam2DEulerBernoulli.get_e_size()

                elif self.elements[0].finite_element_type == Beam2DShearDeformable.get_element_type():
                    self.element_e_n = Beam2DShearDeformable.get_e_size()

                else:
                    print Warning, "Element type not supported! Defined element type: ", self.elements[0].finite_element_type

                self.node_e_n = self.element_e_n / 2

            else:
                print Warning, "Elements list is empty!"

        #   predefine empty array
        e = np.array([])
        self.node_dim = []

        #   only nodal coordinates vector of first node of every element is added to the vector of nodal coordinates of the mesh
        #   and from the last element in a mesh also the absolute nodal coordinates of the last node is also added to the nodal coordinates of the mesh
        for i, element in enumerate(self.elements):
            #   every element except last
            if i != self.n_e - 1:
                #   check slope discontinuity
                #   first node is not checked for node discontinuity
                if i == 0:
                    e_i = element.evaluate_e()[0:element.node_e_n]

                    #   append vector of nodal coordinates to the vector of mesh nodal coordinates
                    e = np.r_[e, e_i]
                    # print "e 1. node =", e
                    self.node_dim.append(len(e_i))

                else:
                    #   no slope discontinuity
                    if element.theta[2] == self.elements[i-1].theta[2]:
                        e_i = element.evaluate_e()[0:element.node_e_n]
                        #   append vector of nodal coordinates to the vector of mesh nodal coordinates
                        e = np.r_[e, e_i]
                        self.node_dim.append(len(e_i))

                    #   slope discontinuity
                    else:
                        # slopde_discontinuity = SlopeDiscontinuity()
                        # if self._parent is not None:
                        #     print "testing =", self._parent

                        e_i = element.evaluate_e()[0:2]
                        #   structural nodal gradients
                        r_X = np.array([1., 0.])
                        r_Y = np.array([0., 1.])
                        e_i = np.append(e_i, r_X)
                        e_i = np.append(e_i, r_Y)

                        #   append vector of nodal coordinates to the vector of mesh nodal coordinates
                        e = np.r_[e, e_i]
                        self.node_dim.append(len(e_i))

            #   only first element if only one elements is in a mesh (this is only executed if only one element is in a mesh)
            elif self.n_e == 1:
                e_i = element.evaluate_e()
                e_i_i = e_i[0:element.node_e_n]
                e_i_j = e_i[element.node_e_n::]

                #   append vector of nodal coordinates to the vector of mesh nodal coordinates
                e = np.r_[e, e_i_i]
                e = np.r_[e, e_i_j]
                self.node_dim.append(len(e_i_i))
                self.node_dim.append(len(e_i_j))

            #   only last element (if only one element is in mesh the following code is not executed)
            else:
                # print "LAST ELEMENT"
                #   check slope discontinuity
                #   no slope discontinuity
                # print "element.theta[2] =", element.theta[2]
                # print "self.elements[i-1].theta[2] =", self.elements[i-1].theta[2]
                if np.isclose(element.theta[2], self.elements[i-1].theta[2], rtol=1e-06, atol=1e-08, equal_nan=False):
                # if element.theta[2] == self.elements[i-1].theta[2]:
                #     print "NO SLOPE DISCONTIUITY"
                    e_i = element.evaluate_e()
                    # print "e_i =", len(e_i)
                    #   append vector of nodal coordinates to the vector of mesh nodal coordinates
                    e = np.r_[e, e_i]
                    # print "e =", len(e)
                    self.node_dim.append(len(e_i[0:element.node_e_n]))
                    self.node_dim.append(len(e_i[element.node_e_n::]))

                #   slope discontinuity
                else:
                    # print "SLOPE DISCONTIUITY"
                    #   predefined empty array
                    e_i = np.array([])
                    #   vector of nodal coordinates of an element
                    _e_i = element.evaluate_e()
                    #   coordinates of the first node of the element
                    e_i_N1 = _e_i[0:2]
                    e_i = np.append(e_i, e_i_N1)
                    #   structural nodal gradients
                    r_X = np.array([1., 0.])
                    r_Y = np.array([0., 1.])
                    e_i = np.append(e_i, r_X)
                    e_i = np.append(e_i, r_Y)
                    #   append vector of nodal coordinates to the vector of mesh nodal coordinates
                    e = np.r_[e, e_i]
                    self.node_dim.append(len(e_i))

                    #   coordinates of the last node of the element
                    e_i = _e_i[element.node_e_n::]

                    #   append vector of nodal coordinates to the vector of mesh nodal coordinates
                    e = np.r_[e, e_i]
                    self.node_dim.append(len(e_i))

        #   flatten numpy array
        self.e = np.array(e).flatten()
        # print "self.e OUT =", len(self.e)
        # print "self.node_dim =", self.node_dim
        #   number of nodal coordinates
        self.n_NC = len(self.e)

        return self.e

    def evaluate_r(self, e, node_id, element_id=None, ksi_element=None):
        """

        :param e:
        :param body_id:
        :return:
        """
        rP_ij = np.zeros(2)
        if element_id is None and ksi_element is None:
            for element in self.elements:
                # print "element =", element
                if node_id in element.node_id_list:
                    # print "node_id inside!!!!"
                    # print "e =", e
                    e_i = element.evaluate_e_i(e_b=e)
                    # e_i = np.array([0.1, 0., 1., 0., 0.15, 0., 1., 0.])
                    for i, _node_id in enumerate(element.node_id_list):
                        # print "i =", i
                        if _node_id == node_id:
                            rP_ij = element.evaluate_r(ksi=i, e_i=e_i)
                            break

        else:
            element = self.elements[element_id]
            e_i = element.evaluate_e_i(e_b=e)
            rP_ij = element.evaluate_r(ksi=ksi_element, e_i=e_i)

        return rP_ij

    def evaluate_dr(self, de, element_id=None, ksi_element=None):
        """

        :param de:
        :param element_id:
        :param ksi_element:
        :return:
        """
        element = self.elements[element_id]
        de_i = element.evaluate_e_i(e_b=de)
        drP = element.evaluate_r(ksi=ksi_element, e_i=de_i)

        return drP

    def e_node(self, e, node_id):
        """

        :return:
        """
        e_node = e[node_id * self.node_e_n:node_id * self.node_e_n + self.node_e_n]
        return e_node

    def evaluate_S(self, e, node_id):
        """
        Function evaluates shape function of
        :param e:           of a mesh
        :param node_id:     node id
        :return:
        """
        S_i = None
        for element in self.elements:
            # print "element =", element
            if node_id in element.node_id_list:
                e_i = element.evaluate_e_i(e_b=e)
                # e_i = np.array([0.1, 0., 1., 0., 0.15, 0., 1., 0.])
                for i, _node_id in enumerate(element.node_id_list):
                    # print "i =", i
                    if _node_id == node_id:
                        S_i = element._evaluate_S(i)
                        # print "rP_ij =", rP_ij
                        return S_i
                        break

        return S_i


    def check_slope_discontinuity(self, parent):
        """

        :return:
        """
        for i, (node_dim_i, node_id) in enumerate(zip(self.node_dim, self.node_id_list)):

            if node_dim_i == 6 and i < len(self.node_id_list):
                slopde_discontinuity = SlopeDiscontinuity(self._parent.body_id, node_id_i=node_id, node_id_j=i, parent=parent)
                self.slope_discontinuities.append(slopde_discontinuity)

    def evaluate_dq0(self):
        """

        :return:
        """
        self.de = np.zeros_like(self.e)
        return self.de

    def e_j(self, node_id, e=None):
        """
        Function returns a vector of nodal coordinate vectors at node j
        :param node_id:
        :return:
        """
        if e is not None:
            self.e = e

        #   nodal coordinates of a node j with node id
        e_j = self.e[sum(self.node_dim[0:node_id]):sum(self.node_dim[0:node_id + 1])]

        return e_j

    def evaluate_ksi_element_id_node_id_at_position(self, L):
        """

        :return:
        """
        ksi = None
        node_id = None
        element_id = None

        self.L0 = 0.
        for element in self.elements:
            self.L0 += element.L

            if self.L0 >= L:
                ksi = (self.L0 - L) / element.L
                element_id = element.element_id
                #   if length is equal to length to last node of element
                if self.L0 == L:
                    node_id = element.node_id_list[-1]

        return ksi, element_id, node_id

    def plot_mesh(self, ax=None, save=False, show=False):
        """

        :return:
        """
        if ax is None:
            fig = plt.figure(figsize=(6, 5),
                             dpi=100,
                             facecolor='w',
                             edgecolor='k')

            ax = plt.subplot(111, aspect="equal")
            ax.ticklabel_format(style='sci', axis='both')

        for element in self.elements:
            element.plot(ax)

        plt.grid(True)
        if show:
            plt.show()

        # plt.xlim([-1, 3])
        # plt.ylim([-1, 3])

        if save:
            plt.savefig(self._name+".png")
            plt.clf()

    def plot_geometry(self, t=None, e=None, e_b=None, ax=None, save=False, show=False, color=np.array([0, 0, 0])):
        """
        Function plot undeformed and deformed geometry
        :return:
        """
        if e is not None:
            self.e = e
        else:
            self.e = self.evaluate_q()

        if ax is None:
            fig = plt.figure(figsize=(6, 5),
                             dpi=100,
                             facecolor='w',
                             edgecolor='k')

            ax = plt.subplot(111, aspect="equal")
            ax.ticklabel_format(style='sci', axis='both')

            dxy = 1.5
            plt.xlim([-dxy, dxy])
            plt.ylim([-dxy, dxy])

            #   grid
            plt.grid(True)

        #   evaluate mesh geometry at current state
        print "e_b =", e_b
        self.evaluate_geometry_nodes(e_b=e_b)

        #   plot geometry nodes
        plt.plot(np.array(self.geometry_nodes)[:, 0], np.array(self.geometry_nodes)[:, 1], color=color)

        #   simulation info
        d = 0.1
        x = np.array(self.geometry_nodes)[-1, 0]
        y = np.array(self.geometry_nodes)[-1, 1]

        # if self.step is not None:
        #     if mpl.rcParams['text.usetex']:
        #         plt.text(x, y, "step="+str(self.step).zfill(4), fontdict=None)
        #     else:
        #         plt.figtext(x, y, "step="+str(self.step).zfill(4), fontdict=None)
        #
        # if t is not None:
        #     if mpl.rcParams['text.usetex']:
        #         plt.text(x, y, r"\textit{t =}", fontdict=None)
        #     else:
        #         plt.text(x, y, "t =", fontdict=None)

        if save:
            filename = "t_" + str(t).zfill(2) + ".png"
            fig.savefig( filename)
            fig.clf()

    def write_K_to_file(self):
        """
        Method writes stiffness matrix to .dat file
        :return:
        """
        K = self.evaluate_K()

        filename = "stiffness_matrix.dat"

        np.savetxt(K, filename)
        print "Stiffness matrix K saved to file %s"%filename

    def evaluate_kinetic_energy(self, de):
        """

        :return:
        """
        self._kinetic_energy = 0.5 * reduce(np.dot, [de.T, self.M, de])
        return self._kinetic_energy

    def evaluate_potential_energy(self, e_b, g):
        """

        :return:
        """
        self._potential_energy = 0.
        for element in self.elements:
            e_i = element.evaluate_e_i(e_b)
            Ep_i = element.evaluate_potential_energy(e_i, g)
            self._potential_energy += Ep_i

        return self._potential_energy

    def evaluate_elastic_strain_energy(self, e):
        """

        :return:
        """
        self.evaluate_K(e)
        self._strain_energy = 0.5 * reduce(np.dot, [e.T, self.K, e])
        return self._strain_energy

if __name__ == "__main__":
    # filename = "mesh_2d_2elements.mesh"
    # filename = "mesh_2d_frame_2elements.mesh"
    # filename = "example_mesh0.mesh"
    # filename = "example_mesh.mesh"
    # filename = "example_mesh_2elements.mesh"#OK
    # filename = "frame_mesh_2d_2elements.mesh"
    filename = "beam_mesh_2d_2elements.mesh"
    mesh = Mesh()
    mesh.scale = 1.
    mesh.element_type = "beam 2D euler-bernoulli"
    #   assign geometry and mass properties
    mesh.cross_section.type = "rectangle"
    mesh.cross_section.B = .1
    mesh.cross_section.H = .5

    mesh.cross_section.evaluate()

    mesh.density = 7800.
    mesh.module_of_elasticity = 2.07E+11

    mesh.read_file(file_abs_path=filename)
    # mesh.evaluate_geometry_nodes()
    print "nodes0-OK=", mesh.nodes0
    # mesh.sort_nodes()
    print "nodes-NOK=", mesh.nodes
    print mesh.evaluate_q()
    # print "e_mesh =", len(mesh.e)
    # print "node_dim =", mesh.node_dim
    # print "e =", mesh.e

    # for i, element in enumerate(mesh.elements):
    #     print i, element, element.node_id_i, element.node_id_j
    #     B = element.evaluate_B()
    #     T = element.evaluate_T()
    #     M = element.evaluate_M()
    #
    #     TB = np.dot(T, B)

    # M = mesh.evaluate_M()
    # K = mesh.evaluate_K(mesh.evaluate_q())
    # print "K ="
    # print K
    # print "mesh.e ="
    # print mesh.e
    # print "Q_S ="
    # print np.dot(K, mesh.e)

    # fig = plt.figure(figsize=(6, 5),
    #                  dpi=100,
    #                  facecolor='w',
    #                  edgecolor='k')
    #
    # ax = plt.subplot(111, aspect="auto")
    # ax.ticklabel_format(style='sci', axis='both')
    # mesh.plot_mesh(ax=ax, show=True)

    # Q_g = mesh.evaluate_Q_g(g=np.array([0., -9.83]))
    # print "Q_g ="
    # print Q_g





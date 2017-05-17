"""
Created on 17. nov. 2013

@author: lskrinjar
"""
import time
import copy
import itertools
import warnings
from pprint import pprint


import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import vtk


from MBD_system.A import A_matrix
from MBD_system.Ai_ui_P import Ai_ui_P_vector
from MBD_system.q2R_i import q2R_i
from MBD_system.q2theta_i import q2theta_i
from MBD_system.transform_cs import gcs2lcs_z_axis
from bounding_box_data_container import DataContainerAABB
from MBD_system.transform_cs import uP_lcs2gcs

np.set_printoptions(precision=4, threshold=None, edgeitems=100, linewidth=1000, suppress=False, nanstr=None, infstr=None)


class AABB2D(object):  # AABBItem or object - for testing
    """
    classdocs
    https://en.wikipedia.org/wiki/K-d_tree
    """
    __level = 1
    __id = itertools.count(1)

    def __init__(self, nodes_LCS=[], normals_LCS=[], tangents_LCS=[], angles_LCS=[], min_LCS=None, max_LCS=None, z_dim=0, parent_body=None, id_sub=1, profile_type="closed", _type=None, parent=None):
        # super(AABBItem, self).__init__(name="AABB" + parent._name + "_body=", parent=parent)
        """
        Constructor of BB (bounding box) from nodes of a geometry
        Args:
            nodes_LCS     - nodes body matrix (columns = 3, rows = N, dtype = numpy array) in local coordinate system of a body
            normals_LCS   - normals body matrix (columns = 3, rows = N, dtype = numpy array)
            min_LCS       - 
            max_LCS       - 
            z_dim         - dimension in z direction (a plane) where can be contact between two bodies
            parent_body   - a pointer to body object that this AABB item is part of AABB tree of a body
            id_sub        - 
            parent        - 
            visibility    - 
        Constructs:
            tree of BB objects
        """
        self._parent = parent
        self.children = []
        self._typeInfo = "aabb"
        self._parent_body = parent_body

        self.id_sub = id_sub
        self.id = self.__id.next()
        self._type = _type
        self.subdivide_direction = None

        #   dimension of data
        self.len_nodes, self.nD = np.shape(nodes_LCS)

        #    nodes in AABB as empty list - predefined empty list
        # self.len_nodes, self.nD = np.shape(nodes_LCS)
        # #   2D
        # if self.nD != 3:
        #     #    add zeros to third column - z component
        #     self.nodes_LCS = np.c_[nodes_LCS, np.zeros([self.len_nodes, 1])]
        #     if normals_LCS != []:
        #         print np.shape(normals_LCS)
        #         self.normals_LCS = np.c_[normals_LCS, np.zeros([self.len_nodes, 1])]
        #     else:
        #         self.normals_LCS = normals_LCS
        #
        # #   3D
        # else:

        #   data type
        self.profile_type = profile_type
        #   nodes, normals, tangents in LCS
        self.nodes_LCS = nodes_LCS
        self.normals_LCS = normals_LCS
        self.tangents_LCS = tangents_LCS
        self.angles_LCS = angles_LCS
        
        self.min_LCS = min_LCS
        self.max_LCS = max_LCS

        #   nodes, normals in AABB
        #   LCS
        self.nodes_in_AABB_LCS = []
        self.normals_in_AABB_LCS = []
        self.tangents_in_AABB_LCS = []

        #   GCS
        self.nodes_in_AABB_GCS = []
        self.normals_in_AABB_GCS = []
        self.tangents_in_AABB_GCS = []

        #   geometry properties
        self.x_center = self.y_center = self.z_center = 0.
        self.x_center_GCS = self.y_center_GCS = self.z_center_GCS = 0.

        #   visualization properties
        self._visible = False
        self._VBO_created = False
        self.vtk_actor = None
        self.color = np.array([1.0, 0.0, 0.0], dtype="float32")

        # print "self._parent._typeInfo =", self._parent._typeInfo
        # print "self._parent_body._typeInfo =", self._parent_body._typeInfo
        #    this if is only for the example with matplotlib display of AABB
        if self._parent==None and self._parent_body==None:
            self.level = self.__level
            self.id_path = [1]
            
            self.body_id = 0
            self.max_nodes_LCS_in_AABB = 3
            self.min_dimensions_of_AABB = .0001
            self.max_level = 6
            self.z_dim = z_dim

            self.color = np.random.rand(3)
            
        #   when building AABBtree from contact - at first step
        elif self._parent._typeInfo == "contact" and self._parent_body._typeInfo == "body":
            self.profile_type = "opened"
            self.level = self.__level
            self.id_path = [1]
            self.body_id = self._parent_body.body_id
            self.min_dimensions_of_AABB = 0.0001
            self.max_level = 6
            self._parent_body = parent_body
            self.max_nodes_LCS_in_AABB = 4

            #    color
            self.color = np.array([1.0, 0.0, 0.0], dtype="float32")
            # self.z_dim = self._parent.z_dim
            # print "************************************************"
            # print self.nodes_LCS
            # print "self._parent_body._name =", self._parent_body._name
            # print "self._parent.z_dim =", self._parent.z_dim
            # print "self._parent_body.R[2] =", self._parent_body.R[2]
            # print "self._parent_body.CM_CAD_LCS[2] =", self._parent_body.CM_CAD_LCS[2]
            # print "self._parent_body.CAD_LCS_GCS[2] =", self._parent_body.CAD_LCS_GCS[2]
            self.z_dim = gcs2lcs_z_axis(+self._parent_body.CM_CAD_LCS[2], self._parent.z_dim)
            # self.z_dim = self._parent.z_dim
            # print "self.z_dim =", self.z_dim
            # np.savetxt(self._parent_body._name+"_LCS_2D.txt", self.nodes_LCS)

        #   when building AABBtree from contact - each step after first in recursion
        elif self._parent._typeInfo == "aabb" and self._parent_body == None:
            self.profile_type = "opened"

            self.level = self._parent.level + 1
            __id_path = copy.copy(self._parent.id_path)
            self.id_path = __id_path
            self.id_path.append(id_sub)
            self.body_id = self._parent.body_id
            self.max_nodes_LCS_in_AABB = self._parent.max_nodes_LCS_in_AABB
            self.min_dimensions_of_AABB = self._parent.min_dimensions_of_AABB
            self.max_level = self._parent.max_level

            self.id_sub = id_sub
            self._divided_parent_list = self._parent.children
            
            self._parent_body = self._parent._parent_body
            self.z_dim = self._parent.z_dim

            #    color
            self.color = np.random.rand(3)
            # self.z_dim = gcs2lcs_z_axis(self._parent_body.CM_CAD_LCS[2], self._parent.z_dim)

        else:
            raise ValueError, "None of the conditions are met!"

        # print "level-id =", self.level,"-",self.id
        #    first AABB root object of each geometry is first constructed and than additional parameters are assigned to this object
        #    if additional parameters are defined in data files they are assigned to this root object at level 1
        #    at level 1 method construct() is called explicitly, when building lower AABBs (higher level) method construct() is called implicitly in __init__() 
        if self.level == 1:
            pass
        else:
            self.construct()

    def set_vtk_data(self, interactor, source):
        """

        :return:
        """
        self.frame = vtk.vtkOutlineFilter()
        self.frame.SetInputData(source)

        #   vtk mapper
        self.vtk_mapper = vtk.vtkPolyDataMapper()
        self.vtk_mapper.SetInputConnection(self.frame.GetOutputPort())

        #   vtk actor
        self.vtk_actor = vtk.vtkActor()
        self.vtk_actor.SetMapper(self.vtk_mapper)
        self.vtk_actor.SetPosition(self._parent_body.R)
        self.vtk_actor.SetOrientation(np.rad2deg(self._parent_body.theta))

    def update_vtk_data(self, q):
        """

        :param q:
        :return:
        """
        self.vtk_actor.SetPosition(self._parent_body.R)
        self.vtk_actor.SetOrientation(np.rad2deg(self._parent_body.theta))

    def add_additional_parameters(self, key, value):
        """
        Args:
            key - parameter name (str)
            value - value of parameter name (num or str)
        """
        setattr(self, key, value)

    def _get_properties(self):
        """

        :return:
        """
        print "---------------------"
        print "body :", self._parent_body._name
        pprint(vars(self))

    def construct(self):
        """
        Function creates AABB tree object
        """
        #    opengl properties
        #    visibility
        self._visible = True
        #    vbo status
        self._VBO_created = False

        if self._parent is None or self._parent_body is None:
            #    R(x, y)
            self.R = np.zeros(3)
            #    theta
            self.theta = np.zeros(3)
        else:
            #    R(x, y)
            self.R = self._parent_body.R  # self._parent_body.R
            #    theta
            self.theta = self._parent_body.theta  # self._parent_body.theta

        # self.R = np.zeros(3)
        # self.theta = np.zeros(3)

        #    constructed AABB boundary min-max points
#         self._boundary(_nodes=self.nodes_LCS, _max=self.max_LCS, _min=self.min_LCS)
        self._boundary(_nodes=self.nodes_LCS, _max=self.max_LCS, _min=self.min_LCS)
            
        #    construct frame geometry (points array, indices array for VBO and IBO)
        self._frame_geometry_LCS = self.frame_geometry_LCS()

        #   frame in GCS
        q_i = np.zeros(3)
        if self._parent is not None:
            if self._parent._parent is not None:
                if hasattr(self._parent._parent._parent, "evaluate_q"):
                    q_i = self._parent._parent._parent.evaluate_q()

        self.frame_geometry_GCS(q_i)

        #   status
        self.constructed = True
        #    this is part of the the recursion to subdivide nodes and create
        if (len(self.nodes_LCS) <= self.max_nodes_LCS_in_AABB) or (self._dimensions < self.min_dimensions_of_AABB).any() or (self.level == self.max_level):# or (self.level == 1):
            self.nodes_in_AABB_LCS = np.array(self.nodes_LCS)
            self.normals_in_AABB_LCS = np.array(self.normals_LCS)
            self.tangents_in_AABB_LCS = np.array(self.tangents_LCS)
            self.angles_in_AABB_LCS = np.array(self.angles_LCS)
            self.phi = []
            # print "----------------------------------------"
            # print str(self.level) + "-" + str(self.id)
            # print "self.nodes_in_AABB_LCS ="
            # print self.nodes_in_AABB_LCS
        else:
            self._subdivide(nodes=self.nodes_LCS, normals=self.normals_LCS, tangents=self.tangents_LCS, angles=self.angles_LCS)

    def _subdivide(self, nodes, normals, tangents, angles):
        """
        Subdivide AABB
        Args:
            nodes - nodes (np.array) to be divided
            normals - nodes (np.array) to be divided
        Returns:
            None
        Raises:
            None
        """
        self.sub_data_list = []
        #    true-false(TF) list of nodes that are sub nodes
        _sub_nodes_min_TF = nodes[:, self.divide_direction] <= self.divide_node[self.divide_direction]
        _sub_nodes_max_TF = nodes[:, self.divide_direction] >= self.divide_node[self.divide_direction]

        #   join to one list to process in for loop
        _sub_nodes_TF_list = [_sub_nodes_min_TF, _sub_nodes_max_TF]

        #    subdivide in x direction
        if self.dx >= self.dy:
            self.subdivide_direction = "x"
            #    first min, then max
            for __id, _sub_nodes_TF in enumerate(_sub_nodes_TF_list):
                _sub_nodes = nodes[_sub_nodes_TF]
                #    sub box boundary
                #    min sub box in x direction
                if (_sub_nodes_TF == _sub_nodes_min_TF).all():
                    __type = "dx_min"
                    # print "min x dir"
                    _x_min = _sub_nodes[:, 0].min()
                    _x_max = np.delete(nodes[_sub_nodes_max_TF][:, 0], self.divide_node).min()
                    # _y_min =
                #    max sub box in x direction
                if (_sub_nodes_TF == _sub_nodes_max_TF).all():
                    __type = "dx_max"
                    # print "max x dir"
                    _x_min = np.delete(nodes[_sub_nodes_min_TF][:, 0], self.divide_node).max()
                    _x_max = _sub_nodes[:, 0].max()

                #    y boundary is common to min and max sub box
                _y_max = self.y_max
                _y_min = nodes[:, 1].min()
                _sub_normals = self._sub_vectors(normals, _sub_nodes_TF)
                _sub_tangents = self._sub_vectors(tangents, _sub_nodes_TF)
                _sub_angles = self._sub_vals(angles, _sub_nodes_TF)
                # print "_sub_nodes", type(_sub_nodes)
                # print _sub_nodes
                data_containter = DataContainerAABB(_id=__id, _nodes=_sub_nodes, _normals=_sub_normals, _tangents=_sub_tangents, _angles=_sub_angles, x_min=_x_min, x_max=_x_max, y_min=_y_min, y_max=_y_max, _type = __type)
                self.sub_data_list.append(data_containter)

        #    subdivide in y direction
        if self.dy > self.dx:
            self.subdivide_direction = "y"
            #    first min, then max
            for __id, _sub_nodes_TF in enumerate(_sub_nodes_TF_list):
                _sub_nodes = nodes[_sub_nodes_TF]
                #    sub box boundary
                #    min sub box in y direction
                if (_sub_nodes_TF == _sub_nodes_min_TF).all():
                    __type = "dy_min"
                    _x_min = _sub_nodes[:, 0].min()
                    _y_min = _sub_nodes[:, 1].min()
                    # _y_min = np.delete(nodes[_sub_nodes_min_TF][:, 1], np.where(nodes[_sub_nodes_min_TF][:, 1] == self.divide_node[1])[0][0]).min()
                    # _y_max = np.delete(nodes[_sub_nodes_max_TF][:, 1], self.divide_node).min()
                    _y_max = _sub_nodes[:, 1].max()
                    # _y_max = np.delete(nodes[_sub_nodes_max_TF][:, 1], np.where(nodes[_sub_nodes_max_TF][:, 1] == self.divide_node[1])[0][0]).min()

                    xy_min = np.array([self.x_min, _y_min])

                #    max sub box in y direction
                if (_sub_nodes_TF == _sub_nodes_max_TF).all():
                    __type = "dy_max"
                    # print "max y dir"
                    _x_min = np.delete(nodes[_sub_nodes_min_TF][:, 0], np.where(nodes[_sub_nodes_min_TF][:, 0] == self.divide_node[0])[0][0]).min()

                    _y_min = np.delete(nodes[_sub_nodes_min_TF][:, 1], np.where(nodes[_sub_nodes_min_TF][:, 1] == self.divide_node[1])[0][0]).max()
                    # _y_min = np.delete(nodes[_sub_nodes_max_TF][:, 1], np.where(nodes[_sub_nodes_max_TF][:, 1] == self.divide_node[1])[0][0]).min()
                    # _y_min = _sub_nodes[:, 1].min()
                    _y_max = _sub_nodes[:, 1].max()
                    # _y_max = np.delete(nodes[_sub_nodes_max_TF][:, 1], np.where(nodes[_sub_nodes_max_TF][:, 1] == self.divide_node[1])[0][0]).min()

                    xy_min = np.array([_x_min, _y_min])

                #   check if node [x_min, y_min] is in nodes and if not add it to max sub box
                if xy_min.tolist() not in _sub_nodes.tolist() and xy_min.tolist() in nodes.tolist():
                    #   get indices
                    _indx = np.where(xy_min == nodes)

                    #   indices by row and column
                    _rows, _cols = _indx

                    #   change bool item at index position
                    _sub_nodes_TF[_rows[0]] = True

                    #   change sub nodes
                    _sub_nodes = nodes[_sub_nodes_TF]

                #    common to min and max sub box
                # _x_max = _sub_nodes[:, 0].max()
                _x_max = self.x_max
                _x_min = self.x_min
                _sub_normals = self._sub_vectors(normals, _sub_nodes_TF)
                _sub_tangents = self._sub_vectors(tangents, _sub_nodes_TF)
                _sub_angles = self._sub_vals(angles, _sub_nodes_TF)
                # print "_y_min = ", _y_min
                data_containter = DataContainerAABB(_id=__id, _nodes=_sub_nodes, _normals=_sub_normals, _tangents=_sub_tangents, _angles=_sub_angles, x_min=_x_min, x_max=_x_max, y_min=_y_min, y_max=_y_max, _type=__type)
                self.sub_data_list.append(data_containter)

        #   check connectivity
        #   when division in y direction
        if self.subdivide_direction == "y":
            #   check if last node is in max sub division
            if nodes[-1, :].tolist() in self.sub_data_list[1]._nodes.tolist():
                pass

            #   check if first node is in min sub division
            if nodes[0, :].tolist() in self.sub_data_list[0]._nodes.tolist():
                #   check if last node in min sub division is equal to the first node of max sub AABB
                if np.equal(self.sub_data_list[1]._nodes[0, :], self.sub_data_list[0]._nodes[-1, :]).all(axis=0).all():
                    print "last node in min sub AABB is equal to first of max sub AABB"
                # if np.equal(nodes[0, :], self.sub_data_list[0]._nodes[-1, :]).all(axis=0).any():
                # if nodes[0, :] == self.sub_data_list[0]._nodes[-1, 0]:
                    pass
                else:
                    print "last node in min sub AABB is not equal to first of max sub AABB", self.id
                    if self.id == 6:
                        print "last node in min sub AABB should be =", self.sub_data_list[1]._nodes[0, :]
                        print "max sub AABB nodes ="
                        print self.sub_data_list[1]._nodes
                    # if self.id == 6:
                    #     print "nodes before reposition"
                    #     print "first node in nodes =", nodes[0, :]
                    #     print "nodes in min sub section"
                    #     print self.sub_data_list[0]._nodes

                    #   first remove it from current position and append it to last position
                    _nodes = np.delete(self.sub_data_list[0]._nodes, nodes[0, :], axis=0)
                    _nodes = np.vstack((_nodes, nodes[0, :]))
                    self.sub_data_list[0]._nodes = _nodes

                    # if self.id == 6:
                    #     print "nodes after reposition"
                    #     print "nodes in min sub section"
                    #     print self.sub_data_list[0]._nodes

                # print "nodes in min sub section (after)"
                # print self.sub_data_list[0]._nodes

                #   check if first node in min sub division is equal to first node in parent AABB
                if np.equal(nodes[0, :], self.sub_data_list[0]._nodes[0, :]).all(axis=0).all():
                    print "first node in min sub AABB is equal to first of parent AABB", "id ="
                    print "first node in min sub AABB =", self.sub_data_list[0]._nodes[0, :]
                    print "first node of parent AABB =", nodes[0, :]
                else:
                    print "first node in min sub AABB is not equal to first of parent AABB", self.id

            #   check if first node in
            if nodes[0, :].tolist() in self.sub_data_list[0]._nodes.tolist():
                pass

        #   when division in x direction
        if self.subdivide_direction == "x":
            #   if node is already in list check the position and reorder if needed
            if nodes[0, :].tolist() in self.sub_data_list[0]._nodes.tolist():
                pass
                # print "node is inside, check for position"
                # print "nodes ="
                # print nodes
                # if all(self.sub_data_list[0]._nodes[-1, :] == ):
                # print "self.sub_data_list[0]._type =", self.sub_data_list[0]._type
                #
                # if all():
                #     self.sub_data_list[0]._nodes.tolist()
            else:
                pass

        #    initiate 2 children
        for sub_data_obj in self.sub_data_list:
            sub_AABB = AABB2D(nodes_LCS=sub_data_obj._nodes, normals_LCS=sub_data_obj._normals, tangents_LCS=sub_data_obj._tangents, angles_LCS=sub_data_obj._angles_LCS, min_LCS=sub_data_obj._min, max_LCS=sub_data_obj._max, id_sub=sub_data_obj._id, parent=self, profile_type="opened", _type=sub_data_obj._type)
            self.children.append(sub_AABB)

    def _check_ratio(self, x_min, x_max, y_min, y_max):
        """
        Function checks if ratio is inside/outside of limits 
        """
        #    dimensions
        x_dim = x_max - x_min
        y_dim = y_max - y_min
    
        #   check width/height ratio
        _ratio = min(x_dim, y_dim) / max(x_dim, y_dim)
        _factor = 0.1

        if _ratio < _factor:
            #    modify x dimension. if x dimension is too small
            if x_dim < _factor * y_dim:
                y_min = y_min
                y_max = y_max
                dx = _factor * y_dim
    
                if x_min == 0:
                    x_min = -0.5 * dx
                else:
                    x_min = x_min - np.sign(x_min) * dx
                    
                if x_max == 0:
                    x_max = +0.5 * dx
                else:
                    x_max = x_max + np.sign(x_max) * dx

            #    modify y dimension, if y dimension is too small
            if y_dim < _factor * x_dim:
                x_min = x_min
                x_max = x_max
                dy = _factor * x_dim
                
                if y_min == 0:
                    y_min = -0.5 * dy
                else:
                    if np.sign(y_min) == 1:
                        y_min = y_min - dy
                    if np.sign(y_min) == -1:
                        y_min = y_min - dy
                    
                if y_max == 0:
                    y_max = +0.5 * dy
                else:
                    if np.sign(y_max) == 1:
                        y_max = y_max + dy
                    if np.sign(y_max) == -1:
                        y_max = y_max + dy

        return x_min, x_max, y_min, y_max

    def _sub_vectors(self, vectors, sub_indices):
        """
        Function returns sub vector based on boolean indices
        """
        if vectors!=[]:
            _sub_vectors = vectors[sub_indices]
        else:
            _sub_vectors = []
            
        return _sub_vectors

    def _sub_vals(self, vals, sub_indices):
        """

        Args:
            vals:
            sub_indices:

        Returns:
        """
        if vals != []:
            _sub_vals = [vals for vals, sub_indx in zip(vals, sub_indices) if sub_indx]
            # _sub_vals = itertools.compress(vals, sub_indices)  # vectors[sub_indices]
        else:
            _sub_vals = []

        return _sub_vals

    def _boundary(self, _nodes, _max=None, _min=None):
        """
        Construct boundary frame
        """
        if _min is None:
            self._min = np.array(_nodes).min(axis=0)
        else:
            self._min = _min

        if _max is None:
            self._max = np.array(_nodes).max(axis=0)
        else:
            self._max = _max

        #    get components
        #    max-min in x,y,z dimensions in LCS
        #    2D
        if len(self._min) == 2:
            [self.x_min, self.y_min] = self._min
            [self.x_max, self.y_max] = self._max
            self.z_min = self.z_dim
            self.z_max = self.z_dim
            
        #    3D    
        if len(self._min) == 3:
            [self.x_min, self.y_min, self.z_min] = self._min
            [self.x_max, self.y_max, self.z_max] = self._max

        self.z_min = self.z_dim
        self.z_max = self.z_dim

        #    check x/y to height dimension ratio
        self.x_min, self.x_max, self.y_min, self.y_max = self._check_ratio(self.x_min, self.x_max, self.y_min, self.y_max)
        #    calculate dimensions of a AABB frame
        self._frame_dimensions(self.x_min, self.x_max, self.y_min, self.y_max)

        #    frame center
        self._center(_nodes)

        #    get divide node - node nearest to the center    
        self._divide_node(_nodes)

    def _frame_dimensions(self, x_min, x_max, y_min, y_max, z_min=0, z_max=0):
        """
        Calculate dimensions of AABB frame and direction in which the subdivision is made
        """
        #    dimensions of a bounding box (frame)
        self.dimensions = np.array([x_max - x_min, y_max - y_min])  # max - _min
        
        #    2D 
        if len(self.dimensions) == 2:
            [self.dx, self.dy] = self.dimensions
            self.dz = 0
        #    3D
        elif len(self.dimensions) == 3:
            [self.dx, self.dy, self.dz] = self.dimensions
        
        #    remove zero elements from dimension to get 2D data
        self._dimensions = self.dimensions[np.nonzero(self.dimensions)]

        #    get max and min dimensions of a frame
        max_dim = self.dimensions.max()
        min_dim = self.dimensions.min()

        #    divide direction is index of max element of dimensions array
        self.divide_direction = np.where(max_dim == self.dimensions)[0][0]

        #    get position of zero element in array
        zero_elem_pos = np.where(self.dimensions == 0)[0]
        #    if a zero element is found it is replaced with a numerical small value
        if len(zero_elem_pos) != 0:
            self.dimensions[zero_elem_pos] = 1E-100

    def _center(self, _nodes):
        """
        Construct center of AABB
        """
        #   center node is computed as geometry center of all points
        self.center = np.sum(_nodes, axis=0) / len(_nodes)
        #   get components (x,y,z)
        self.x_center, self.y_center = self.center

        #   coordinates in GCS
        [self.x_center_GCS, self.y_center_GCS] = self.R[0:2] + Ai_ui_P_vector(self.center, self.theta[2])#[self.x_center, self.y_center, self.z_center]

    def _divide_node(self, _array):
        """
        Function calculates (selects) divide node, that is the node that is nearest
        to the center
        """
        self.divide_node = _array[(np.abs(_array[:, self.divide_direction] - self.center[self.divide_direction])).argmin()]

    def _contact_geometry_GCS(self, q):
        """
        Transform nodes from LCS to GCS
        :param q:   vector of MBD system
        """
        #   nodes in GCS
        self.nodes_in_AABB_GCS = q2R_i(q, self.body_id) + Ai_ui_P_vector(self.nodes_in_AABB_LCS, q2theta_i(q, self.body_id))
        #   normals in GCS
        self.normals_in_AABB_GCS = Ai_ui_P_vector(self.normals_in_AABB_LCS, q2theta_i(q, self.body_id))
        #   tangents in GCS
        self.tangents_in_AABB_GCS = Ai_ui_P_vector(self.tangents_in_AABB_LCS, q2theta_i(q, self.body_id))
        # print "self.normals_in_AABB_GCS =", self.normals_in_AABB_GCS[0,:]
        # print "self.tangents_in_AABB_GCS =", self.tangents_in_AABB_GCS[0,:]

    def _frame_GCS(self, q=np.zeros(3)):
        """

        :param q:
        :return:
        """

        [self.x_min_GCS, self.y_min_GCS, self.z_min_GCS] = np.append(q[0:2], self.z_dim) + Ai_ui_P_vector(np.array([self.x_min, self.y_min, self.z_min]), q[2])
        [self.x_max_GCS, self.y_max_GCS, self.z_max_GCS] = np.append(q[0:2], self.z_dim) + Ai_ui_P_vector(np.array([self.x_max, self.y_max, self.z_max]), q[2])

        print '_frame_GCS='
        print [self.x_min_GCS, self.y_min_GCS, self.z_min_GCS]
        print [self.x_max_GCS, self.y_max_GCS, self.z_max_GCS]

    def frame_geometry_LCS(self):
        """
        Construct frame geometry points.
        """
        #    construct if contact data is 2D or 3D
        if self.nD == 2:
            frame_nodes_LCS = np.array([
                                         [self.x_min, self.y_min, self.z_min],  # T1
                                         [self.x_max, self.y_min, self.z_min],  # T2
                                         [self.x_max, self.y_max, self.z_min],  # T3
                                         [self.x_min, self.y_max, self.z_min],  # T4
                                         ], dtype='float32')
        else:
            frame_nodes_LCS = np.array([
                                         [self.x_min, self.y_min, self.z_min],  # T1
                                         [self.x_max, self.y_min, self.z_min],  # T2
                                         [self.x_max, self.y_max, self.z_min],  # T3
                                         [self.x_min, self.y_max, self.z_min],  # T4
                                         [self.x_min, self.y_min, self.z_max],  # T5
                                         [self.x_max, self.y_min, self.z_max],  # T6
                                         [self.x_max, self.y_max, self.z_max],  # T7
                                         [self.x_min, self.y_max, self.z_max]  # T8
                                         ], dtype='float32')
        

        self.VBO_data = np.array(frame_nodes_LCS.flatten(), dtype='float32')

        if self.level == 1:
            if self.nD == 2:
                self.frame_indx = np.array([0, 1, 2, 3, 0], dtype='int32')
            else:
                self.frame_indx = np.array([0, 1, 2, 3,
                                    4, 5, 6, 7,
                                    0, 4, 1, 5, 2, 6, 3, 7], dtype='int32')
        else:
            self.frame_indx = self._parent.frame_indx

        return frame_nodes_LCS

    def has_children(self):
        """
        Function returns:
        true - if object has children
        false- if object does not have children (if object does not have children it has nodes)
        """
        if self.children == []:
            return False
        else:
            return True

    def has_nodes(self):
        """
        Function returns:
        true - if AABB has nodes inside
        false - if AABB does not have nodes inside
        """
        if len(self.nodes_LCS) > 0:
            return True
        else:
            return False

    def has_nodes_in_AABB(self):
        """
        Function returns:
        true - if AABB has nodes inside
        false - if AABB does not have nodes inside
        """
        if len(self.nodes_in_AABB_LCS) > 0:
            return True
        else:
            return False

    def number_of_nodes_inside(self):
        """
        Function returns number of nodes inside AABB boundary.
        """
        return len(self.nodes_in_AABB_LCS)

    def frame_geometry_GCS(self, q_i):
        """
        Function updates frame geometry to use it with contact detection - AABB overlap on CPU.
        :param q:   body vector
        """
        R_i = q_i[0:2]
        theta_i = q_i[2]
        self._frame_geometry_GCS = R_i + Ai_ui_P_vector(self._frame_geometry_LCS[:, 0:2], theta_i)

        self.x_min_GCS, self.y_min_GCS = np.amin(self._frame_geometry_GCS, axis=0)
        self.x_max_GCS, self.y_max_GCS = np.amax(self._frame_geometry_GCS, axis=0)

        self.update_frame_geometry_center_2D(q_i)

    def update_frame_geometry(self, q=None):
        """

        :param q:
        :return:
        """
        self._frame_GCS(q)

    def update_frame_geometry_center_2D(self, q_i):
        """
        Function updates positio vector to center of AABB frame
        """
        R_i = q_i[0:2]
        theta_i = q_i[2]
        [self.x_center_GCS, self.y_center_GCS] = R_i + Ai_ui_P_vector(np.array([self.x_center, self.y_center]), theta_i)

    def update_nodes_GCS_in_AABB_2D(self, q):
        """
        Function updates nodes GCS only when two AABB overlap. New values of attribute R 
        are already calculated in function update_frame_geometry(q), because this function 
        executes after function update_frame_geometry(q) and after condition of AABB overlap is achieved.
        """
        xy_nodes_GCS_in_AABB = np.zeros([len(self.nodes_in_AABB_LCS), 2])
        #    calculate xy coordinates of nodes inside AABB
        R = q2R_i(q, self.body_id)
        theta = q2theta_i(q, self.body_id)
        for i, node_LCS in enumerate(self.nodes_in_AABB_LCS):
            xy_nodes_GCS_in_AABB[i, :] = R + Ai_ui_P_vector(node_LCS[0:2], theta)

        self.nodes_in_AABB_GCS = xy_nodes_GCS_in_AABB

    def transform_LCS2GCS(self, q, uP):
        """
        Function transforms LCS data in GCS data
        :return:
        """
        #    calculate xy coordinates of nodes inside AABB
        R = q2R_i(q, self.body_id)
        theta = q2theta_i(q, self.body_id)

        rP = []
        for uPi in uP:
            rPi = R + Ai_ui_P_vector(uPi[0:2], theta)
            rP.append(rPi)

        return rP

    def _update_color(self, color):
        self.color = color

    def print_nodes_LCS(self):
        print "Nodes in LCS"
        print self.nodes_in_AABB_LCS

    def print_nodes_GCS(self):
        if self.nodes_in_AABB_GCS==[]:
            q = self._parent._parent._parent.evaluate_q()
            self.update_nodes_GCS_in_AABB_2D(q)

        print "Nodes in GCS"
        print self.nodes_in_AABB_GCS

    def plot_2D(self, _ax=None, color=None):
        """
        
        """
        if color is None:
            color = self.color

        for child in self.children:
            child.plot_2D(_ax=_ax)

        if _ax is None:
            fig = plt.figure(figsize=(10, 8),
                             dpi=300,
                             facecolor='w',
                             edgecolor='k')
            ax = plt.subplot(111, aspect="auto")
            ax.ticklabel_format(style='sci', axis='both')
        else:
            ax = _ax

        #   plot frame
        if self.has_nodes_in_AABB():
            _fill = True
        else:
            _fill = False

        #   evaluate GCS of AABB
        if self._parent_body is not None:
            q_i = self._parent_body.get_q()
            self.frame_geometry_GCS(q_i=q_i)

        frame = matplotlib.patches.Rectangle(np.array([self.x_min_GCS, self.y_min_GCS]),
                                             self.x_max_GCS - self.x_min_GCS,
                                             self.y_max_GCS - self.y_min_GCS,
                                             angle=0.0,
                                             fill=_fill,
                                             alpha=0.2,
                                             color=color)
        ax.add_patch(frame)

        #   plot center of AABB
        ax.plot(self.x_center_GCS, self.y_center_GCS, color=color, marker="x")
        if self._type is None:
            self._type = "_"
        #   plot id of AABB at center
        if self.subdivide_direction is None:
            ax.text(self.x_center_GCS * 1.1, self.y_center_GCS * 1.1, "ID: " + str(self.id) + "(" + self._type + ")", color="black")
        else:
            ax.text(self.x_center_GCS * 1.1, self.y_center_GCS * 1.1, "ID: " + str(self.id) + " division in: " + self.subdivide_direction.upper(), color="black")

        #   plot children id
        for i, child in enumerate(self.children):
            ax.text(self.x_center_GCS * 1.1, self.y_center_GCS * 1.1 + ((i + 1)), "ID_" + str(i) + ": " + str(child.id), color=child.color)

        #   plot AABBs with and without nodes inside
        q = np.zeros(3)
        if self._parent_body is not None:
            q = self._parent_body._parent._parent.evaluate_q()
        if self.has_nodes_in_AABB():
            self._contact_geometry_GCS(q)
            _nodes = self.nodes_in_AABB_GCS

        else:
            _nodes = np.array(self.transform_LCS2GCS(q, self.nodes_LCS))

        #   for visualization if profile is closed, first node is appended to array to connect last node with first node when plotting node data
        if self.profile_type == "closed":
            _nodes = np.vstack((_nodes, _nodes[0, :]))

        #   plot nodes and lines
        plt.plot(_nodes[:, 0],
                 _nodes[:, 1],
                 color=color,
                 linestyle="-",
                 marker="s")

        #   plot coordinates
        for i, node in enumerate(_nodes):
            plt.text(node[0] * 1.2, node[1] * 1.2, "ID = " + str(i), fontsize=10)
            #   coordinates as text
            plt.text(node[0] * 1.4, node[1] * 1.4, "[" + str(round(node[0], 3)) + "," + str(round(node[1], 3)) + "]", fontsize=8)

        #   plot info
        plt.text(0.1, 0.9,'Coordinate vectors are in mm.',
                 horizontalalignment='left',
                 verticalalignment='center',
                 transform = ax.transAxes,
                 fontsize=8)

        #   save figure to file
        if _ax is None:
            # plt.xlim([1.2 * _nodes.min(), 1.2 * _nodes.max()])
            # plt.ylim([1.2 * _nodes.min(), 1.2 * _nodes.max()])

            # xlim = [-12E-3, +12E-3]
            # ylim = [-12E-3, +12E-3]
            #
            # plt.xlim(xlim)
            # plt.ylim(ylim)

            plt.grid(True)
            filename = "plot_AABB_id_" + str(self.id).zfill(4) + ".png"
            if __name__ == "__main__":
                pass
            else:
                print "Plot saved to file %s" % filename
            plt.savefig(filename)
            # fig.show()
            fig.clf()

    def plot_3D(self, ax, _color):
        """
        Function plots the AABB boundary if AABB has nodes inside.
        """
        
        if self.has_nodes():
            T1 = np.array([self.x_min, self.y_min, self.z_min])
            T2 = np.array([self.x_max, self.y_min, self.z_min])
            T3 = np.array([self.x_max, self.y_max, self.z_min])
            T4 = np.array([self.x_min, self.y_max, self.z_min])
            T5 = np.array([self.x_min, self.y_min, self.z_max])
            T6 = np.array([self.x_max, self.y_min, self.z_max])
            T7 = np.array([self.x_max, self.y_max, self.z_max])
            T8 = np.array([self.x_min, self.y_max, self.z_max])
    
            self.matrix = np.array([
                   T1,
                   T2,
                   T3,
                   T4,
                   T1,
                   T5,
                   T6,
                   T2,
                   T6,
                   T7,
                   T3,
                   T7,
                   T8,
                   T4,
                   T8,
                   T5])
            _color = np.random.rand(3)
#             color_ = np.array([1, 0, 0])
            
            ax.plot_wireframe(self.matrix[:, 0], self.matrix[:, 1], self.matrix[:, 2], rstride=10, cstride=10, color=_color)
            
            ax.scatter(self.x_center, self.y_center, self.z_center, color=_color, marker="x")
            
            ax.text(self.x_center, self.y_center, self.z_center, str(self.level) + "-" + str(self.id), color='black')

    def get_tree_items(self, tree_item):
        """
        
        """
        if hasattr(tree_item, "R"):
            print tree_item.R
        else:
            print tree_item.id, tree_item.id_path
        print "*******************************"
        
        for tree_item_ in tree_item.children:
            self.get_tree_items(tree_item_)

    def __del__(self):
        """

        """

    def _show(self):
        """

        """
        if self._visible:
            self._visible = False
        else:
            self._visible = True


def remove_duplicate_nodes(body):
    """
    Function removes duplicate nodes
    """
    nodes = body.geom.geom_data.vertices
    normals = body.geom.geom_data.normals
    
    nodes_ = np.ascontiguousarray(nodes).view(np.dtype((np.void, nodes.dtype.itemsize * nodes.shape[1])))
    none_, idx = np.unique(nodes_, return_index=True)
    
    nodes_ = nodes[np.sort(idx)]
    normals_ = np.repeat(normals, repeats=3, axis=0)[np.sort(idx)]
    return nodes_, normals_


class AABBlist(object):
    def __init__(self, filename):
        self.nodes = np.array(np.loadtxt(filename, delimiter=','))  # nodes_cube_problem, nodes_cube, nodes_cube_problem_nok, nodes_cube_problem_ok, nodes_cube_problem_ok
        self.tree = AABB2D(nodes_LCS=self.nodes)
        self.tree.max_level = 4
        self.tree.min_dimensions_of_AABB = 10
        self.tree.construct()

        self._list = []
        
    def get_children(self, obj):
        self._list.append(obj)
        if obj.has_children():
            for child in obj.children:
                self.get_children(child)
        

if __name__ == "__main__":
    warnings.filterwarnings("ignore")
    # filename = "latch_lever_LCS_2D.txt", "circle.dat"
    # np.savetxt("circle.dat", data, fmt='%.4e', delimiter='\t')

    # filename = "circle.dat"
    # filename = "circle_n=10.txt"
    filename = "circle_n_20.txt"
    nodes = np.loadtxt(filename, delimiter=",")
    # print "nodes ="
    # print nodes
    AABB = AABB2D(nodes_LCS=nodes, profile_type="closed")
    AABB.max_level = 4
    #   OK to level 3, level 4 NOK
    AABB.max_nodes_LCS_in_AABB = 3
    AABB.construct()

    AABB.plot_2D()
    print "Finished succesfully!"



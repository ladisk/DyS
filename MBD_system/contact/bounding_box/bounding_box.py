"""
Created on 17. nov. 2013

@author: lskrinjar
"""
import warnings
import time
import copy
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from pprint import pprint
import ctypes
from OpenGL import *
from OpenGL.GL import *
from OpenGL.GL.shaders import *
from OpenGL.GLU import *
from dxfwrite import DXFEngine as dxf
import itertools


from bounding_box_data_container import DataContainerAABB
from MBD_system.A import A_matrix
from MBD_system.MBD_system_items import AABBItem
from MBD_system.transform_cs import gcs2lcs_z_axis
from MBD_system.q2R_i import q2R_i
from MBD_system.q2theta_i import q2theta_i
from MBD_system.Ai_ui_P import Ai_ui_P_vector

np.set_printoptions(precision=4, threshold=None, edgeitems=100, linewidth=1000, suppress=False, nanstr=None, infstr=None)


class AABB(object):  # AABBItem or object - for testing
    """
    classdocs
    """
    __level = 1
    __id = itertools.count(1)

    def __init__(self, nodes_LCS=[], normals_LCS=[], min_LCS=None, max_LCS=None, z_dim=None, parent_body=None, id_sub=1, parent=None, visibility=False, _type=None):
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

        #    nodes in AABB as empty list - predefined empty list
        self.len_nodes, self.nD = np.shape(nodes_LCS)

        if self.nD != 3:
            #    add zeros to third column - z component
            self.nodes_LCS = np.c_[nodes_LCS, np.zeros([self.len_nodes, 1])]
            if normals_LCS != []:
                self.normals_LCS = np.c_[normals_LCS, np.zeros([self.len_nodes, 1])]
            else:
                self.normals_LCS = normals_LCS
        else:
            self.nodes_LCS = nodes_LCS
            self.normals_LCS = normals_LCS

        self.min_LCS = min_LCS
        self.max_LCS = max_LCS

        self.nodes_LCS_in_AABB = []        
        self.normals_LCS_in_AABB = []
        self.nodes_GCS_in_AABB = []
        self.normals_GCS_in_AABB = []

        self._visible = False
        self._VBO_created = False

        #    this if is only for the example with matplotlib display of AABB
        if self._parent == None and self._parent_body == None:
            self.level = self.__level
            self.id_path = [1]
            
            self.body_id = 0
            self.max_nodes_LCS_in_AABB = 3
            self.min_dimensions_of_AABB = .001
            self.max_level = 6
            self.z_dim = 0
            
        #   when building AABBtree from contact - at first step
        elif self._parent._typeInfo == "contact" and self._parent_body._typeInfo == "body":
            self.level = self.__level
            self.id_path = [1]
            self.body_id = self._parent_body.body_id
            self.min_dimensions_of_AABB = None
            self.max_level = None
            self._parent_body = parent_body
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
            np.savetxt(self._parent_body._name+"_LCS_2D.txt", self.nodes_LCS)


        #   when building AABBtree from contact - each step after first in recursion
        elif self._parent._typeInfo == "aabb" and self._parent_body == None:
            self.level = self._parent.level + 1
            __id_path = copy.copy(self._parent.id_path)
            self.id_path = __id_path
            self.id_path.append(id_sub)
            
            self.max_nodes_LCS_in_AABB = self._parent.max_nodes_LCS_in_AABB
            self.min_dimensions_of_AABB = self._parent.min_dimensions_of_AABB
            self.max_level = self._parent.max_level
            
            self._parent_body = self._parent
            self.id_sub = id_sub
            self._divided_parent_list = self._parent.children
            
            self._parent_body = self._parent._parent_body
            self.z_dim = self._parent.z_dim
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
        Function creates AABB object
        """
        # print "construct()"
        #    opengl properties
        #    AABB frame color
        self.color_GL = np.array([1.0, 0.0, 0.0], dtype="float32")
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
        self._frame_geometry()

        #   status
        self.constructed = True
        #    this is part of the the recursion to subdivide nodes and create
        if (len(self.nodes_LCS) <= self.max_nodes_LCS_in_AABB) or (self._dimensions < self.min_dimensions_of_AABB).any() or (self.level == self.max_level) or (self.level == 1):
            self.nodes_LCS_in_AABB = np.array(self.nodes_LCS)
            self.normals_LCS_in_AABB = np.array(self.normals_LCS)
            self.phi = []
        else:
            self._subdivide(nodes=self.nodes_LCS, normals=self.normals_LCS)

    def _subdivide(self, nodes, normals):
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
        #    true-false(TF) list of nodes that are min sub nodes
        _sub_nodes_min_TF = nodes[:, self.divide_direction] <= self.divide_node[self.divide_direction]
        _sub_nodes_max_TF = nodes[:, self.divide_direction] >= self.divide_node[self.divide_direction]

        #   join to one list to process in for loop
        _sub_nodes_TF_list = [_sub_nodes_min_TF, _sub_nodes_max_TF]


        #    subdivide in x direction
        if self.dx >= self.dy:
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
                    # _y_min = self._parent.y_min

                #    y boundary is common to min and max sub box
                _y_max = self.y_max
                _y_min = nodes[:, 1].min()
                _sub_normals = self.__sub_normals(normals, _sub_nodes_TF)
                # print "_sub_nodes", type(_sub_nodes)
                # print _sub_nodes
                data_containter = DataContainerAABB(_id=__id, _nodes=_sub_nodes, _normals=_sub_normals, x_min=_x_min, x_max=_x_max, y_min=_y_min, y_max=_y_max, _type = __type)
                self.sub_data_list.append(data_containter)
        
        
        #    subdivide in y direction
        if self.dy > self.dx:
            #    first min, then max
            for __id, _sub_nodes_TF in enumerate(_sub_nodes_TF_list):
                _sub_nodes = nodes[_sub_nodes_TF]
                #    sub box boundary
                #    min sub box in y direction
                if (_sub_nodes_TF == _sub_nodes_min_TF).all():
                    __type = "dy_min"
                    _x_min = _sub_nodes[:, 0].min()
                    # _y_min = _sub_nodes[:, 1].min()
                    _y_min = np.delete(nodes[_sub_nodes_min_TF][:, 1], np.where(nodes[_sub_nodes_min_TF][:, 1] == self.divide_node[1])[0][0]).min()
                    _y_max = np.delete(nodes[_sub_nodes_max_TF][:, 1], self.divide_node).min()
                    _y_max = np.delete(nodes[_sub_nodes_max_TF][:, 1], np.where(nodes[_sub_nodes_max_TF][:, 1] == self.divide_node[1])[0][0]).min()

                #    max sub box in y direction
                if (_sub_nodes_TF == _sub_nodes_max_TF).all():
                    __type = "dy_max"
                    # print "max y dir"
                    _x_min = np.delete(nodes[_sub_nodes_min_TF][:, 0], np.where(nodes[_sub_nodes_min_TF][:, 0] == self.divide_node[0])[0][0]).min()

                    _y_min = np.delete(nodes[_sub_nodes_min_TF][:, 1], np.where(nodes[_sub_nodes_min_TF][:, 1] == self.divide_node[1])[0][0]).max()
                    # _y_min = np.delete(nodes[_sub_nodes_max_TF][:, 1], np.where(nodes[_sub_nodes_max_TF][:, 1] == self.divide_node[1])[0][0]).min()

                    _y_max = _sub_nodes[:, 1].max()
                #    common to min and max sub box
                # _x_max = _sub_nodes[:, 0].max()
                _x_max = self.x_max
                _x_min = self.x_min
                _sub_normals = self.__sub_normals(normals, _sub_nodes_TF)
                # print "_y_min = ", _y_min
                data_containter = DataContainerAABB(_id=__id, _nodes=_sub_nodes, _normals=_sub_normals, x_min=_x_min, x_max=_x_max, y_min=_y_min, y_max=_y_max, _type = __type)
                self.sub_data_list.append(data_containter)

        #    initiate 2 children
        for sub_data_obj in self.sub_data_list:
            sub_AABB = AABBTree(nodes_LCS=sub_data_obj._nodes, normals_LCS=sub_data_obj._normals, min_LCS=sub_data_obj._min, max_LCS=sub_data_obj._max, id_sub=sub_data_obj._id, parent=self, _type = sub_data_obj._type)
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
                _y_min = y_min
                _y_max = y_max
                dx = _factor * y_dim
    
                if x_min == 0:
                    _x_min = -0.5 * dx
                else:
                    _x_min = x_min - np.sign(x_min) * dx
                    
                if x_max == 0:
                    _x_max = +0.5 * dx
                else:
                    _x_max = x_max + np.sign(x_max) * dx

            #    modify y dimension, if y dimension is too small
            if y_dim < _factor * x_dim:
                _x_min = x_min
                _x_max = x_max
                dy = _factor * x_dim
                
                if y_min == 0:
                    _y_min = -0.5 * dy
                else:
                    if np.sign(y_min) == 1:
                        _y_min = y_min - dy
                    if np.sign(y_min) == -1:
                        _y_min = y_min - dy
                    
                if y_max == 0:
                    _y_max = +0.5 * dy
                else:
                    if np.sign(y_max) == 1:
                        _y_max = y_max + dy
                    if np.sign(y_max) == -1:
                        _y_max = y_max + dy
        return _x_min, _x_max, _y_min, _y_max


    def __sub_normals(self, normals, sub_nodes_TF):
        """
        
        """
        if normals != []:
            _sub_normals = normals[sub_nodes_TF]
        else:
            _sub_normals = []
            
        return _sub_normals

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
        self.x_center, self.y_center, self.z_center = self.center

        #   coordinates in GCS
        [self.x_center_GCS, self.y_center_GCS, self.z_center_GCS] = self.R + Ai_ui_P_vector(self.center, self.theta[2])#[self.x_center, self.y_center, self.z_center]

    def _divide_node(self, _array):
        """
        Function calculates (selects) divide node, that is the node that is nearest
        to the center
        """
        self.divide_node = _array[(np.abs(_array[:, self.divide_direction] - self.center[self.divide_direction])).argmin()]

    def _lcs2gcs(self):
        """
        Transform LCS to GCS
        """
        [self.x_min_GCS, self.y_min_GCS, self.z_min_GCS] = [self.x_min, self.y_min, self.z_min] + self.R  # add theta
        [self.x_max_GCS, self.y_max_GCS, self.z_max_GCS] = [self.x_max, self.y_max, self.z_max] + self.R  # add theta

    def _frame_geometry(self):
        """
        Construct frame geometry points.
        """
        #    construct if contact data is 2D or 3D
        if self.nD == 2:
            self.frame_nodes_LCS = np.array([
                                         [self.x_min, self.y_min, self.z_min],  # T1
                                         [self.x_max, self.y_min, self.z_min],  # T2
                                         [self.x_max, self.y_max, self.z_min],  # T3
                                         [self.x_min, self.y_max, self.z_min],  # T4
                                         ], dtype='float32')
        else:
            self.frame_nodes_LCS = np.array([
                                         [self.x_min, self.y_min, self.z_min],  # T1
                                         [self.x_max, self.y_min, self.z_min],  # T2
                                         [self.x_max, self.y_max, self.z_min],  # T3
                                         [self.x_min, self.y_max, self.z_min],  # T4
                                         [self.x_min, self.y_min, self.z_max],  # T5
                                         [self.x_max, self.y_min, self.z_max],  # T6
                                         [self.x_max, self.y_max, self.z_max],  # T7
                                         [self.x_min, self.y_max, self.z_max]  # T8
                                         ], dtype='float32')
        

        self.VBO_data = np.array(self.frame_nodes_LCS.flatten(), dtype='float32')


        if self.level == 1:
            if self.nD == 2:
                self.frame_indx = np.array([0, 1, 2, 3, 0], dtype='int32')
            else:
                self.frame_indx = np.array([0, 1, 2, 3,
                                    4, 5, 6, 7,
                                    0, 4, 1, 5, 2, 6, 3, 7], dtype='int32')
        else:
            self.frame_indx = self._parent.frame_indx

        self._lcs2gcs()

    def _create_VBO(self):
        """
        Create VBO - vertex buffer object
        """
        #    generate a new VBO and get the associated with vbo_id
        self.vbo_id = GLuint()
        self.vbo_id = glGenBuffers(1)


        #    bind name to buffer
        glBindBuffer(GL_ARRAY_BUFFER, self.vbo_id)
        #    test if is buffer
        if glIsBuffer(self.vbo_id) == 1:
            self._VBO_created = True
        else:
            self._VBO_created = False
            raise Error, "VBO not created!"

        #    add VBO data to buffer
        data_size_in_bytes = arrays.ArrayDatatype.arrayByteCount(self.VBO_data)
        glBufferData(GL_ARRAY_BUFFER, data_size_in_bytes, self.VBO_data, GL_STATIC_DRAW)
        glBindBuffer(GL_ARRAY_BUFFER, 0)

    def _create_IBO(self):
        """
        Create IBO - index buffer object
        """
        #    generate a new IBO and get the associated with ibo_id
        self.ibo_id = GLuint()
        self.ibo_id = glGenBuffers(1)
        
        #    bind name to buffer
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, self.ibo_id)
        #    test if is buffer
        if glIsBuffer(self.ibo_id) == 1:
            self.VBO_created = True
        else:
            self.VBO_created = False
            raise Error, "VBO not created!"

        #    generate a buffer for the indices for ibo_id
        data_size_in_bytes = arrays.ArrayDatatype.arrayByteCount(self.frame_indx)
        glBufferData(GL_ELEMENT_ARRAY_BUFFER, data_size_in_bytes, self.frame_indx, GL_STATIC_DRAW)
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0)

    def create_VBO(self):
        """
        
        """
        if self.level == 1:
            #    only one IBO - index buffer object for an AABB tree
            self._create_IBO()
        self._create_VBO()

    def paintGL_VBO(self):
        """
        
        """
        if self.level == 1:
            pass
        else:
#            if hasattr(self._parent, "_typeInfo"):
#                if self._parent._typeInfo == self._typeInfo:
            self.ibo_id = self._parent.ibo_id
        
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, self.ibo_id)
        
        
        glBindBuffer(GL_ARRAY_BUFFER, self.vbo_id)
        
        attrib_pos = 0
        glEnableVertexAttribArray(attrib_pos)

        glVertexAttribPointer(attrib_pos,  # attribute
                                3,  # number of elements per vertex, here (x,y,z)
                                GL_FLOAT,  # the type of each element
                                False,  # take our values as-is
                                0,  # no extra data between each position
                                None)  # offset of first element


        glLineWidth(1.)
        glDisable(GL_LIGHTING)

        glDrawElements(GL_LINE_LOOP, 4, GL_UNSIGNED_INT, ctypes.c_void_p(0))
        glDrawElements(GL_LINE_LOOP, 4, GL_UNSIGNED_INT, ctypes.c_void_p(16))
        glDrawElements(GL_LINES, 8, GL_UNSIGNED_INT, ctypes.c_void_p(32))

        glBegin(GL_POINTS)
        self._paintGL_nodes()
        glEnd()

        glEnable(GL_LIGHTING)
        glLineWidth(1.)
        glPointSize(1.0)
        
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0)

        
        glBindBuffer(GL_ARRAY_BUFFER, 0)
        glDisableVertexAttribArray(attrib_pos)


    def _paintGL_nodes(self):
        """
        
        :return:
        """
        # for node in self.nodes_LCS_in_AABB:
        #     glVertex3f(np.float32(node[0]), np.float32(node[1]), np.float32(node[2]))


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
        if len(self.nodes_LCS_in_AABB) > 0:
            return True
        else:
            return False

    
    def number_of_nodes_inside(self):
        """
        Function returns number of nodes inside AABB boundary.
        """
        return len(self.nodes_LCS_in_AABB)


    def update_frame_geometry(self, q_i):
        """
        Function updates frame geometry to use it with contact detection - AABB overlap on CPU.
        q_i - body q vector q = [Rx, Ry, theta]
        """
        # print "update_frame_geometry()"
        R = q_i[0:2]
        theta = q_i[2]

        A_matrix_ = A_matrix(theta)
         
        u_P_min = np.array([self.x_min, self.y_min])
        xy_min = np.dot(A_matrix_, u_P_min)
         
        u_P_max = np.array([self.x_max, self.y_max])
        xy_max = np.dot(A_matrix_, u_P_max)

        [self.x_min_GCS, self.y_min_GCS, self.z_min_GCS] = self._parent_body.R + np.append(xy_min, self.z_min)  # [self.x_min, self.y_min, self.z_min]  # 
        [self.x_max_GCS, self.y_max_GCS, self.z_max_GCS] = self._parent_body.R + np.append(xy_max, self.z_max)  # [self.x_max, self.y_max, self.z_max]  #  
    
        self.update_frame_geometry_center_2D()

    def update_frame_geometry_center_2D(self):
        """
        Function updates positio vector to center of AABB frame
        """
        [self.x_center_GCS, self.y_center_GCS] = [self.x_center, self.y_center] + self.R[0:2] # add theta
        
    def update_nodes_and_normals_GCS_in_AABB_2D(self, q):
        """
        Function updates nodes and normals in AABB
        """
        #   body R, theta
        R_i = q2R_i(q, self._parent_body.body_id)
        theta_i = q2theta_i(q, self._parent_body.body_id)

        self.update_nodes_GCS_in_AABB_2D(R_i, theta_i)
        self.update_normals_GCS_in_AABB_2D(theta_i)
        
        
    def update_nodes_GCS_in_AABB_2D(self, R_i, theta_i):
        """
        Function updates nodes GCS only when two AABB overlap. New values of attribute R 
        are already calculated in function update_frame_geometry(q), because this function 
        executes after function update_frame_geometry(q) and after condition of AABB overlap is achieved.
        """
        xy_nodes_GCS_in_AABB = np.zeros([len(self.nodes_LCS_in_AABB), 2])
        #    calculate xy coordinates of nodes inside AABB
        for i, node_LCS in enumerate(self.nodes_LCS_in_AABB):
            xy_nodes_GCS_in_AABB[i, :] = R_i + Ai_ui_P_vector(node_LCS[0:2], theta_i)

        #    assign z coordinates of nodes inside AABB as planar movement
        z_nodes_GCS_in_AABB = np.array([self.nodes_LCS_in_AABB[:, 2]]).T
        
        self.nodes_GCS_in_AABB = np.hstack((xy_nodes_GCS_in_AABB, z_nodes_GCS_in_AABB))

        
    def update_normals_GCS_in_AABB_2D(self, theta_i):
        """
        
        """
        self.normals_GCS_in_AABB = self.normals_LCS_in_AABB
        for i, normal in enumerate(self.normals_LCS_in_AABB):
            self.normals_GCS_in_AABB[i, :2] = Ai_ui_P_vector(normal[0:2], theta_i)

        
    def plot_2D(self, ax, _color):
        """
        
        """
        T1 = np.array([self.x_min_GCS, self.y_min_GCS])
        T2 = np.array([self.x_max_GCS, self.y_min_GCS])
        T3 = np.array([self.x_max_GCS, self.y_max_GCS])
        T4 = np.array([self.x_min_GCS, self.y_max_GCS])
        
        self.matrix = np.array([
               T1,
               T2,
               T3,
               T4,
               T1])

        # ax.plot(self.matrix[:, 0], self.matrix[:, 1], color=_color)

        # ax.plot(self.x_center_GCS, self.y_center_GCS, color=_color, marker="x")

        frame = matplotlib.patches.Rectangle(T1, self.x_max_GCS-self.x_min_GCS, self.y_max_GCS-self.y_min_GCS, angle=0.0, fill=True, alpha=0.2, color=_color)
        ax.add_patch(frame)
        # _level_id = str(self.level) + "-" + str(self.id)
        # ax.text(self.x_center_GCS, self.y_center_GCS, _level_id, color='black')
        #
        # ax.text(0.01, 0.018, "Type: "+str(self._type), color='black')
        # ax.text(0.01, 0.016, "Nodes inside: "+str(self.has_nodes_in_AABB()), color='black')
        # ax.text(0.01, 0.014, "Number of nodes: "+str(self.number_of_nodes_inside()), color='black')
        # ax.text(0.01, 0.012, "Level-ID: "+_level_id, color='black')
        #
        #
        # if self._parent is None:
        #     _level_id_parent = "None"
        # else:
        #     _level_id_parent = str(self._parent.level) + "-" + str(self._parent.id)
        #
        # ax.text(0.01, 0.010, "Parent (Levep-ID): "+ _level_id_parent, color='black')
        #
        # ax.text(0.01, 0.008, "x (min, max): %8.6f, %8.6f"%(self.x_min_GCS, self.x_max_GCS))
        # ax.text(0.01, 0.006, "y (min, max): %8.6f, %8.6f"%(self.y_min_GCS, self.y_max_GCS))
        # ax.text(0.01, 0.004, "Center (x, y): %8.6f, %8.6f"%(self.x_center_GCS, self.y_center_GCS))
        
        # ax.text(T1[0], T1[1], "min-min", color='black')
        # ax.text(T2[0], T2[1], "max-min", color='black')
        # ax.text(T3[0], T3[1], "max-max", color='black')
        # ax.text(T4[0], T4[1], "min-max", color='black')


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


    def create_VBO_tree(self):
        """
        
        """
        
        self.create_VBO()
            
        for child in self.children:
            child.create_VBO_tree()

    
    def paintGL_VBO_tree(self):
        """
        
        """
        self.paintGL_VBO()
            
        for child in self.children:
            if child._parent_body._name == "ground_":
                pass
            else:
                child.paintGL_VBO_tree()


    def _show_AABB(self):
        """

        """
        if self._visible:
            self._visible = False
        else:
            self._visible = True

    
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
    
    
    def generate_surface_mesh_3D(self):
        """
        Generates surface mesh from triangles based on points of AABB and body geometry data points from file
        """
        None
        
    
    def __del__(self):
        """
        Delete opengl object VBO from VRAM
        """
        if glGetError() == GL_NO_ERROR:
            glDeleteBuffers(1, int(self.vbo_id))
            glDeleteBuffers(1, int(self.ibo_id))


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
        self.nodes_ = np.loadtxt(filename, delimiter='\t')  # nodes_cube_problem, nodes_cube, nodes_cube_problem_nok, nodes_cube_problem_ok, nodes_cube_problem_ok
        self.tree = AABB(nodes_LCS=self.nodes_)
        self.tree.max_level = 4
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

    filename = "circle.dat"

    a = AABBlist(filename)
    a.get_children(a.tree)
    _list = a._list

    fig = plt.figure(num=0, figsize=(10, 6), dpi=160, facecolor='w', edgecolor='k')
    ax = fig.add_subplot(111)
    ax.xlim=[-0.01, +0.02]
#     fig.xlim()"xlim"=[-0.01, +0.02]
    plt.hold(True)


    ax.plot(a.nodes_[:, 0], a.nodes_[:, 1], color="black", linestyle = "-", linewidth=.1)
    for i, AABB in enumerate(_list):
        if True:
            # plt.cla()
            _color = np.random.rand(3)
            # ax.scatter(a.nodes_[:, 0], a.nodes_[:, 1], color="black", marker="o", s=.1)
            # ax.plot(a.nodes_[:, 0], a.nodes_[:, 1], color="black", linestyle = "-")
            # if AABB.has_nodes_in_AABB():
            #     frame = matplotlib.patches.Rectangle((-200,-100), 400, 200, color='yellow')
                # plt.plot(AABB.nodes_LCS_in_AABB[:, 0], AABB.nodes_LCS_in_AABB[:, 1], marker="x", color=_color, ms=10)
                # plt.plot(AABB.nodes_LCS_in_AABB[:, 0], AABB.nodes_LCS_in_AABB[:, 1], "--", color=_color)
            
            ax.set_aspect('equal')
#             plt.xlim(-10, 10)
#             plt.ylim(-10, 10)
            AABB.plot_2D(ax, _color)
#             plt.savefi
#   ("AABB" + str(AABB.level) + "-" + str(AABB.id) + ".png")
#         print "--------------------------"
#         print "level-id =", AABB.level, "-", AABB.id
#         if AABB.level > 1:
#             print "parent.level-id =", AABB._parent.level, "-", AABB._parent.id
#             
#         print "dx =", AABB.dx
#         print "dy =", AABB.dy
#         print "AABB.x_min =", AABB.x_min
#         print "AABB.x_max =", AABB.x_max
#         print "AABB.y_min =", AABB.y_min
#         print "AABB.y_max =", AABB.y_max
#         print "AABB.dim_ratio =", AABB.dim_ratio
#         plt.savefig("plot_" + str(AABB.level)+"-"+str(AABB.id)+".png")
        # plt.pause(0.5)
    plt.savefig("plot_.png")
    print "plot saved!"
#     im_ani = animation.ArtistAnimation(fig, ims, interval=50, repeat_delay=3000, blit=True)
    plt.show()

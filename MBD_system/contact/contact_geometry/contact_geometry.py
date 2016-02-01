"""
Created on 19. feb. 2015

@author: lskrinjar
"""
import os
from pprint import pprint
import numpy as np
from OpenGL import *
from OpenGL.GL import *
from OpenGL.GL.ARB.vertex_buffer_object import *
from OpenGL.GL.shaders import *
from OpenGL.GLU import *

from MBD_system.body.geometry.line.line import Line
from MBD_system.transform_cs import gcs2lcs_z_axis
from MBD_system.transform_cs import cad2cm_lcs


class ContactGeometry(object):
    """
    classdocs
    """

    def __init__(self, _type, body, area=None, filename=None, uP1=None, uP2=None, parent=None):
        """
        Constructor of contact geometry object
        This object has data to be used when performing evaluation if contact between two geometries is present
        The type of geometry is divided into three types:
        (sub)area - if type is area than only a subarea part of entire body geometry is processed and nodes and normals are stored for contact calculation
        edge - if type ie edge than only two points on a body have to be defined an extra option is to simulate a surface roughness of an edge
        body - if type is body than the entire body geometry (data) is used to calculate contact data

        all of this is valid at predefined z dimension (value) as we only work with 2D planar dynamics
        _type - area, edge, body
                
        """
        #    parent is contact object
        self._parent = parent

        #   type
        self._type = _type

        #   pointer to body
        self.body = body

        #   predefined empty list of
        #   nodes
        self.contact_nodes = []
        #   normals
        self.contact_normals = []
        #   edge objects
        self.contact_profile = []

        #   position of contact (plane) in z dimension
        self.z_dim = self._parent.z_dim
        #   in body coordinates
        self.z_dim_lcs = gcs2lcs_z_axis(body.R[2], self.z_dim)

        #   max penetration depth
        self.max_penetration_depth = self.body.max_penetration_depth

        if self._type == "area":
            self._area_contact_profile(area)
        if self._type == "edge":
            self._edge_contact_profile(uP1, uP2)
        if self._type == "body":
            self._body_contour_contact_profile()

        #     #    create 2D profile
        #     self.construct_contact_profile_2D(self.z_dim_lcs, self._min, self._max)
        #
        # if filename is not None and uP1 is not None and uP2 is not None:
        #     self._load_contact_profile(filename)
        #
        #     #    calculate direction of edge vector
        #     self._edge = uP2 - uP1

        # opengl properties
        self.VBO_created = False
        self._visible = True

    def _area_contact_profile(self, _area):
        """

        :return:
        """
        #   area of contact in body CAD CS
        [x_min, x_max, y_min, y_max] = _area

        #   contact area limits
        self._min = [x_min, y_min]
        self._max = [x_max, y_max]

        self.construct_contact_profile_2D(self.z_dim_lcs, _min=self._min, _max=self._max)


    def set_area(self):
        """

        :return:
        """

    def set_edge(self):
        """

        :return:
        """

    def get_nodes(self):
        """

        :return:
        """
        print self.contact_nodes


    def construct_contact_profile_2D(self, z_dim, _min=None, _max=None):
        """
        Function constructs a 2D profile where plane at Z coordinate z_dim intersects mesh of triangles
        :param z_dim:
        :param _min:
        :param _max:
        :return:
        """
        for triangle in self.body.geom.geom_data.triangles:
            if triangle.plane_intersects_triangle(z_dim):

                if _min is not None and _max is not None:
                    if triangle.edge_in_area(_min, _max):
                        for _E in triangle._edge:
                            self.contact_nodes.append(_E)
                            self.contact_normals.append(triangle.normal)

                        # create 2D edge object
                        _edge = Line(triangle._edge[0], triangle._edge[1], triangle.normal)
                        self.contact_profile.append(_edge)

                else:
                    self.contact_nodes.append(triangle._edge[0])
                    self.contact_nodes.append(triangle._edge[1])

                    _edge = Line(triangle._edge[0], triangle._edge[1], triangle.normal)
                    self.profile.append(_edge)

        # remove duplicate nodes
        self.contact_nodes, self.contact_normals = self._remove_duplicate_data(self.contact_nodes, self.contact_normals)

        #    transform nodes from CAD to CM coordinate system
        self.contact_nodes = cad2cm_lcs(self.contact_nodes, self.body.CM_CAD_LCS)

    def get_data(self):
        """
        
        """
        return self.contact_nodes, self.contact_normals

    def _remove_duplicate_data(self, nodes, normals):
        """
        Function removes duplicate nodes
        """
        #   nodes
        _nodes = np.array(nodes)
        #   normals
        _normals = np.array(normals)

        __nodes = np.ascontiguousarray(_nodes).view(np.dtype((np.void, _nodes.dtype.itemsize * _nodes.shape[1])))
        _, idx = np.unique(__nodes, return_index=True)

        _nodes = _nodes[np.sort(idx)]
        _normals = _normals[np.sort(idx)]
        return _nodes, _normals


    def _save_contact_profile(self, nodes, normals=None):
        """
        Save contact profile data, that was created from body geometry.
        """
        if normals is None:
            data = nodes
        else:
            data = np.hstack((nodes, normals))

        np.savetxt(self.body._name + ".prfl", data, "%6.4f", "\t")


    def _load_contact_profile(self, filename):
        """
        Function loads profile data from specified filename
        """
        #    load data with numpy loadtxt
        data = np.array(np.loadtxt(filename, delimiter='\t'))

        #    get size - check if normals are included
        #    else calculate normals
        n, cols = np.shape(data)

        #    if data has only three columns (x,y,z) normals have to be calculated
        if cols == 3:
            self.contact_nodes = data
            self.contact_normals = []

        elif cols == 2:
            pass
        # normals are written in contact profile file (x,y,z,nx,ny,nz) or (x,y,nx,ny)
        else:
            self.contact_nodes = data[:, 0:3]
            self.contact_normals = data[:, 3:]


    def create_VBO(self):
        """

        :return:
        """
        if not self.VBO_created:
            #    generate a new VBO and get the associated vbo_id
            num_of_VBOs = 1

            #    create buffer name
            self.vbo_id = GLuint()
            self.vbo_id = glGenBuffers(num_of_VBOs)

            #    bind name to buffer
            glBindBuffer(GL_ARRAY_BUFFER, self.vbo_id)
            self.N_vertices = len(self.contact_nodes)

            self.vbo_data = np.hstack(
                (self.contact_nodes, np.tile(np.append(self.body.color_GL, 1), (len(self.contact_nodes), 1))))
            self.vbo_data = np.array(self.vbo_data, dtype='float32').flatten()

            #    VBO_data size in bytes
            self.VBO_data_size_in_bytes = arrays.ArrayDatatype.arrayByteCount(self.vbo_data)

            #    add VBO data to buffer
            glBufferData(GL_ARRAY_BUFFER, self.VBO_data_size_in_bytes, self.vbo_data, GL_DYNAMIC_DRAW)
            glBindBuffer(GL_ARRAY_BUFFER, 0)

            self.VBO_created = True

    def paintGL_VBO(self):
        """

        :return:
        """
        if self._visible and self.VBO_created:
            glBindBuffer(GL_ARRAY_BUFFER, self.vbo_id)

            #   pointers
            v_pointer = ctypes.c_void_p(0)
            c_pointer = ctypes.c_void_p(12)

            #    stride in bits (1 float = 4 bits)
            stride_in_bits = 28

            glVertexPointer(3, GL_FLOAT, stride_in_bits, v_pointer)
            glColorPointer(4, GL_FLOAT, stride_in_bits, c_pointer)
            glNormalPointer(GL_FLOAT, 0, ctypes.c_void_p(0))

            glDisable(GL_LIGHTING)
            glDrawArrays(GL_LINE_STRIP, 0, self.N_vertices)
            glDrawArrays(GL_POINTS, 0, self.N_vertices)

            glEnable(GL_LIGHTING)
            glBindBuffer(GL_ARRAY_BUFFER, 0)

"""
Created on 16. nov. 2013

@author: lskrinjar (email: skrinjar.luka@gmail.com)
"""
import itertools
import os
from pprint import pprint


import numpy as np


class ModelLoader(object):
    """
    classdocs
    """
    def __init__(self, filename):
        """
        loads a 3D file type .stl or .obj
        in:
            file_ - file_ = file_name + file_extension
            file_name - 
            file_extension - 
            
        """
        #   filename
        self.filename = filename

        #   node coordinate tolerance
        self.eps = 1E-12

        #   get name adn extension
        (name, self._extension) = os.path.splitext(self.filename)

        #   geometry attributes
        self.vertices = None
        self.normals = None
        self.tangents = None
        self.angles = None

        #    opens the file for reading
        self.geom_file = open(self.filename, "r")

        if self._extension.lower() == '.stl':
            self.load_stl()
        elif self._extension == '.obj':
            self.load_obj()
        elif self._extension == ".txt":
            self.load_txt()
        else:
            raise IOError, "file_extension not correct, check file_extension."

        #   close file
        self.geom_file.close()

    def load_txt(self):
        """
        Function loads data from .txt file of geometry (points) 2D
        :return:
        """
        #   2D
        vertices = np.loadtxt(self.filename, dtype="float32", comments='#', delimiter=",")
        #   add z dimension
        self.vertices = np.hstack((vertices, np.zeros([len(vertices), 1]))) * 1E-3

    def load_stl(self):
        #    predefine parameters
        self.vertices = []
        self.normals = []
        __id_triangle = itertools.count(0)
        
        self.nodes = []
        __id_node = itertools.count(0)
       
        
        for line in self.geom_file:
            strip_line = line.strip()
            values = strip_line.split()
            
            if strip_line.startswith("facet normal"):
                normal = map(float, values[2:5])
                self.normals.append(normal)
                
            elif strip_line.startswith("vertex"):
                vertex = map(float, values[1:4])
                vertex = vertex
                self.vertices.append(vertex)

        #    reshape
        self.vertices = np.array(self.vertices) * 1E-3
        self.vertices[np.abs(self.vertices) < self.eps] = 0
        self.normals = np.array(self.normals)

        #    check if number of vertices is equal to the number of normals
        N_vertices = len(self.vertices)
        N_normals = len(self.normals)

        if N_vertices == 3 * N_normals:
            pass
        else:
            raise IOError, "number of normals and vertices do not a match, check file."

    def get_vertices_3D(self):
        return self.vertices

    def get_vertices_2D(self):
        return self.vertices[:, 0:2] 

    def get_normals_3D(self):
        if np.shape(self.vertices) == np.shape(self.normals):
            return self.normals
        else:
            return np.repeat(self.normals, repeats=3, axis=0)

    def get_normals_2D(self):
        """

        :return:
        """
        if np.shape(self.vertices) == np.shape(self.normals):
            return self.normals[:, 0:2]
        else:
            if self._extension.lower() == '.stl':
                return np.repeat(self.normals, repeats=3, axis=0)[:, 0:2]
            if self._extension.lower() == '.txt':
                return self.normals[:, 0:2]

    def get_tangents_2D(self):
        """

        :return:
        """
        return self.tangents[:, 0:2]

    def get_angles_2D(self):
        """

        :return:
        """
        if self.angles is None:
            self.angles = self.evaluate_angles_2D()

        return self.angles

    def evaluate_angles_2D(self):
        """

        """
        angles = []
        for i, t in enumerate(self.tangents):
            if i==0:
                dotTAN = np.dot(self.tangents[i],-self.tangents[-1])
                angle = np.rad2deg(np.arccos(dotTAN / (np.linalg.norm(self.tangents[i], ord=2) * np.linalg.norm(-self.tangents[-1], ord=2))))

            else:
                dotTAN = np.dot(self.tangents[i-1], -self.tangents[i])
                angle = np.rad2deg(np.arccos(
                dotTAN / (np.linalg.norm(self.tangents[i-1], ord=2) * np.linalg.norm(-self.tangents[i], ord=2))))

            angles.append(angle)

        return angles

    def remove_joined_duplicate(self):
        """
        Function removes duplicate joined nodes and normals
        removes duplicate rows in matrix - equal coordinates and normal 
        """
        #    joined vertices and normals to one matrix to remove only nodes that have equal#
        #    coordinate and normal
        data = np.c_[self.get_vertices_3D(), self.get_normals_3D()]
        __unique = np.unique(data.view([('', data.dtype)] * data.shape[1]))
        unique_data = __unique.view(data.dtype).reshape((__unique.shape[0], data.shape[1]))
        __vertices = unique_data[:, 0:3]
        __normals = unique_data[:, 3:6]
        return __vertices, __normals

    def remove_duplicate_2D(self):
        """
        Remove duplicate nodes that have same x,y coordinates and different z coordinate
        """
        data = np.ascontiguousarray(self.get_vertices_2D())
        __unique, idx = np.unique(data.view([('', data.dtype)] * data.shape[1]), return_index=True)

        _vertices = self.get_vertices_3D()[np.sort(idx)]
        _normals = self.get_normals_3D()[np.sort(idx)]
        return _vertices, _normals

    def remove_duplicate_3D(self):
        """
        Remove duplicate nodes that have same x,y,z coordinates
        :return:
        """
        data = np.ascontiguousarray(self.get_vertices_3D())
        __unique, idx = np.unique(data.view([('', data.dtype)] * data.shape[1]), return_index=True)

        _vertices = self.get_vertices_3D()[np.sort(idx)]
        _normals = self.get_normals_3D()[np.sort(idx)]
        return _vertices, _normals

    def load_obj(self):
        #    predefine parameters
        self.vertices = []
        self.normals = []
        self.texcoords = []
        self.faces = []
        
        material = None
        for line in self.geom_file:
            if line.startswith('#'): continue
            values = line.split()
            if not values: continue
            
            if line.startswith("v "):
                vertex = map(float, values[1:4])
                self.vertices.append(vertex)
                
            elif line.startswith("vn"):
                normal = map(float, values[1:4])
                self.normals.append(normal)
                
            elif line.startswith("vt"):
                texcoord = values[1:3]
                self.texcoords.append(texcoord)
            elif line.startswith("usemtl"):
                material = values[1]
                
        #    reshape
        self.vertices = np.array(self.vertices)
        self.normals = np.array(self.normals)
        
        #    check if number of vertices is equal to the number of normals
        N_vertices = len(self.vertices)
        N_normals = len(self.normals)

        if N_vertices == N_normals:
            pass
        else:
            raise IOError, "number of normals and vertices do not a match, check file."

    def offset_to_CM(self, CM):
        """

        :return:
        """
        self.vertices = self.vertices - CM

    def _create_VBO_array(self, file_extension, colors, GL_primitive_type="triangle", interleaved=True):
        """

        :return:
        """
        #    expand normals matrix, geometry data defines normal per surface-triangle that is defined by 3 points
        if GL_primitive_type == "triangle":
            if file_extension == ".stl":
                if len(self.vertices) != len(self.normals):
                    normals = np.repeat(self.normals, repeats=3, axis=0)

            elif file_extension == ".obj":
                None
            else:
                None

        elif GL_primitive_type == "quad":
            normals = np.repeat(self.normals, repeats=4, axis=0)
            #todo

        elif GL_primitive_type == "line":
            normals = None

        elif GL_primitive_type == "lines":
            normals = None

        else:
            raise ValueError, "GL_primitive_type not correct, check string type"

        #    check if array type is interleaved or not
        if interleaved and GL_primitive_type == "triangle":
            VBO_array = np.hstack((self.vertices, normals, colors)).flatten()

        elif interleaved and GL_primitive_type == "line":
            VBO_array = np.hstack((self.vertices, colors)).flatten()

        elif interleaved and GL_primitive_type == "lines":
            VBO_array = np.hstack((self.vertices, colors)).flatten()

        else:
            VBO_array = np.vstack((self.vertices, normals, colors))

        #   convert data to 32bit float to load to GPU vram
        VBO_array = np.array(VBO_array, dtype='float32')
        return VBO_array

    def construct_contact_profile_2D(self, z_dim, _min=None, _max=None):
        """
        Function constructs a 2D profile where plane at Z coordinate z_dim intersects mesh of triangles
        :param z_dim:
        :param _min:
        :param _max:
        :return:
        """
        self.profile = []
        self.profile_nodes = []
        for triangle in self.triangles:
            if triangle.plane_intersects_triangle(z_dim):

                # if _min is not None and _max is not None:
                #     triangle.edge_in_area(_min, _max)
                #     # self.profile_nodes.append(triangle._p[0])
                #     # self.profile_nodes.append(triangle._p[1])
                #     #
                #     # _edge = Line(triangle._p[0], triangle._p[1], triangle.normal)
                #     # self.profile.append(_edge)
                #
                # else:
                self.profile_nodes.append(triangle._edge[0])
                self.profile_nodes.append(triangle._edge[1])

                _line = Line(triangle._edge[0], triangle._edge[1], triangle.normal)
                self.profile.append(_line)
                
                # for point in triangle._p:
                #     #   if point is not already in list append it
                #     if any((point == node).all() for node in self.profile_nodes):
                #         pass
                #         # self.profile_nodes.append(point)
                #         # self.profile_normals.append(triangle.normal)
                #     else:
                #         self.profile_nodes.append(point)
                #         self.profile_normals.append(triangle.normal)

        # return np.array(self.profile_nodes) * 1E-0, np.array(self.profile_normals)



if __name__ == "__main__":
    # filename = "test_geom.stl"#latch_lever.stl, actuating_lever.stl, test_geom.stl
    filename = "geom.txt"
    a = ModelLoader(filename)
    pprint(vars(a))
    z_dim = 180E-3 #20, 55, 70, 90, 150, 180
    # print "Number of nodes =", 3*len(a.triangles)

    # a.construct_contact_profile_2D(z_dim)
    # print "number of nodes =", len(a.vertices)
    # print "number of nodes in profile ="
    # # print "number of edges =", len(a.profile)
    # print "number of profile nodes", len(a.profile_nodes)
    #
    #
    #
    # fig = plt.figure(num=0, figsize=(8, 6), dpi=120, facecolor='w', edgecolor='k')
    # ax = fig.add_subplot(111)
    # plt.xlim([-.12, 1.2*a.vertices[:, 0].max()])
    # plt.ylim([-.12, 1.2*a.vertices[:, 1].max()])
    # ax.set_aspect('equal')
    # #     fig.xlim()"xlim"=[-0.01, +0.02]
    # plt.hold(True)
    # data = []
    # for edge in a.profile:
    #     _color = np.random.rand(3)
    #     plt.plot(np.array(edge.nodes)[:, 0], np.array(edge.nodes)[:, 1], "-", color=_color, linewidth = 2)
    #
    # #    plot nodes
    # for node in a.vertices:
    #     plt.plot(node[0], node[1], marker="o", markersize = 0.1, color ="black")
    # #    plot lines
    # plt.plot(np.array(a.vertices[:, 0]), np.array(a.vertices[:, 1]), linewidth = 0.05, color =np.array([.8, .8, .8]))
    # #    plot profile nodes
    # for node_ in a.profile_nodes:
    #     plt.plot(node_[0], node_[1], marker="x", markersize = 10, color ="red")
    #
    # print np.array(a.profile_nodes[0]), np.array(a.profile_nodes[1])
        # print edge.nodes
    # plt.plot(np.array(a.profile_nodes)[:,0], np.array(a.profile_nodes)[:,1], "--", color="b")
#     for node in a._profile_nodes:
#         # print node[0], node[1]
#         plt.plot(node[0], node[1], marker="o", color="b", markersize = 1)
#     plt.show()
#
#     name, filetype = filename.split(".")
#     print name
    # print np.array(data)
    # np.savetxt(name+".prfl", data)

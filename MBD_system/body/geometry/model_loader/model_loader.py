"""
Created on 16. nov. 2013

@author: lskrinjar (email: skrinjar.luka@gmail.com)
"""
import os
import itertools
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.animation as animation


from MBD_system.body.geometry.triangle.triangle import Triangle
from MBD_system.body.geometry.node.node import Node
from MBD_system.body.geometry.edge.edge import Edge


class ModelLoader(object):
    """
    classdocs
    """
    def __init__(self, file_):
        """
        loads a 3D file type .stl or .obj
        in:
            file_ - file_ = file_name + file_extension
            file_name - 
            file_extension - 
            
        """
        self.file = file_
        
        self.eps = 1E-12
        
        (filename, file_extension) = os.path.splitext(file_)

        #    opens the file for reading
        self.geom_file = open(self.file, "r")

        if file_extension.lower() == '.stl':
            self.load_stl()
#             print "file_extension: ", file_extension
        elif file_extension == '.obj':
            self.load_obj()
#             print "file_extension:", file_extension
        else:
            raise IOError, "file_extension not correct, check file_extension."

        #   close file
        self.geom_file.close
        
    def load_stl(self):
        #    predefine parameters
        self.vertices = []
        self.normals = []
        self.triangles = []
        __id_triangle = itertools.count(0)
        
        self.nodes = []
        __id_node = itertools.count(0)
       
        
        for line in self.geom_file:
            strip_line = line.strip()
            values = strip_line.split()
            
            if strip_line.startswith("facet normal"):
                normal = map(float, values[2:5])
                self.normals.append(normal)
                #    create triangle object
                _triangle = Triangle(_normal=normal, _id=__id_triangle.next())
                self.triangles.append(_triangle)
                
            elif strip_line.startswith("vertex"):
                vertex = map(float, values[1:4])
                vertex = vertex
                self.vertices.append(vertex)
                
                #    create node object and append it to nodes list
                _node = Node(id_node=__id_node.next(), node=vertex, normal=normal, id_triangle=_triangle._id, parent=_triangle)
                _triangle.add_node(_node)
                
                self.nodes.append(_node)
                

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
        if np.shape(self.vertices) == np.shape(self.normals):
            return self.normals[:, 0:2]
        else:
            return np.repeat(self.normals, repeats=3, axis=0)[:, 0:2]
        
        
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
                #     # _edge = Edge(triangle._p[0], triangle._p[1], triangle.normal)
                #     # self.profile.append(_edge)
                #
                # else:
                self.profile_nodes.append(triangle._edge[0])
                self.profile_nodes.append(triangle._edge[1])

                _edge = Edge(triangle._edge[0], triangle._edge[1], triangle.normal)
                self.profile.append(_edge)
                
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
    filename = "test_geom.stl"#latch_lever.stl, actuating_lever.stl, test_geom.stl
    a = ModelLoader(filename)

    z_dim = 180E-3 #20, 55, 70, 90, 150, 180
    print "Number of nodes =", 3*len(a.triangles)

    a.construct_contact_profile_2D(z_dim)
    print "number of nodes =", len(a.vertices)
    print "number of nodes in profile ="
    # print "number of edges =", len(a.profile)
    print "number of profile nodes", len(a.profile_nodes)
    


    fig = plt.figure(num=0, figsize=(8, 6), dpi=120, facecolor='w', edgecolor='k')
    ax = fig.add_subplot(111)
    plt.xlim([-.12, 1.2*a.vertices[:, 0].max()])
    plt.ylim([-.12, 1.2*a.vertices[:, 1].max()])
    ax.set_aspect('equal')
    #     fig.xlim()"xlim"=[-0.01, +0.02]
    plt.hold(True)
    data = []
    for edge in a.profile:
        _color = np.random.rand(3)
        plt.plot(np.array(edge.nodes)[:, 0], np.array(edge.nodes)[:, 1], "-", color=_color, linewidth = 2)
     

    #    plot nodes
    for node in a.vertices:
        plt.plot(node[0], node[1], marker="o", markersize = 0.1, color ="black")
    #    plot lines
    plt.plot(np.array(a.vertices[:, 0]), np.array(a.vertices[:, 1]), linewidth = 0.05, color =np.array([.8, .8, .8]))
    #    plot profile nodes
    for node_ in a.profile_nodes:
        plt.plot(node_[0], node_[1], marker="x", markersize = 10, color ="red")
        
        
    
    # print np.array(a.profile_nodes[0]), np.array(a.profile_nodes[1])
        # print edge.nodes
    # plt.plot(np.array(a.profile_nodes)[:,0], np.array(a.profile_nodes)[:,1], "--", color="b")
#     for node in a._profile_nodes:
#         # print node[0], node[1]
#         plt.plot(node[0], node[1], marker="o", color="b", markersize = 1)
    plt.show()

    name, filetype = filename.split(".")
    print name
    # print np.array(data)
    # np.savetxt(name+".prfl", data)

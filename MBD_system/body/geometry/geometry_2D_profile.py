# coding=utf-8
"""

created by: lskrinjar
date of creation: 27/03/2016
time of creation: 17:34
"""
import numpy as np
import scipy as sp
import vtk
import logging
import logging.handlers


from MBD_system.body.geometry.geometry_2D import Geometry2D
from MBD_system.t2n import t2n


class ContactGeometry2DProfile(Geometry2D):
    """
    classdocs
    """
    def __init__(self, filename, profile, profile_type="opened", color=np.array([0, 1, 0], dtype="float32"), _dict={}, body=None, parent=None):
        """
        Constructor
        :param parent:
        :return:
        """
        super(ContactGeometry2DProfile, self).__init__(filename, parent=parent)
        #    parent is contact object
        self._parent = parent

        #   pointer to body
        self.body = body

        #   pointer to profile object
        self.profile = profile

        #   visualization properties
        self.vtk_actor = None
        self.color = color

        #   profile type
        #   options: opened or closed (default)
        self.profile_type = profile_type

        #    check if normals and tangets are calculated
        if self.geom_data.tangents is None and self.geom_data.normals is None:
            self._create_normals_and_tangents()

        # if self.geom_data.normals is None:
        #     self._create_normals()

        # print self.geom_data.tangents

    def _create_normals_and_tangents(self):
        """

        :return:
        """
        vertices = self.geom_data.get_vertices_2D()

        #   number of nodes
        n = len(self.geom_data.vertices)

        #   predefine zero matrices
        normals = np.zeros([n, 2])
        tangents = np.zeros([n, 2])
        for i in xrange(0, n-1):
            #   egde vector
            edge = vertices[i+1, :] - vertices[i, :]
            #   unit tangent vector in direction of edge vector
            tangents[i,:] = edge / np.linalg.norm(edge, ord=2)

            #   calculate unit normal of edge
            normals[i, :] = t2n(tangents[i,:])

        #   calculate tangent normal for edge (line) between last and first vertex as they are connected
        if self.profile_type == "closed":
            edge = vertices[0, :] - vertices[-1, :]
            tangents[-1,:] = edge / np.linalg.norm(edge, ord=2)
            normals[-1, :] = t2n(tangents[-1,:])

        self.geom_data.normals = normals
        self.geom_data.tangents = tangents


class ContactGeometry2DClosedProfile(ContactGeometry2DProfile):
    """
    classdocs
    """
    def __init__(self, filename, profile, color=np.array([0, 1, 0], dtype="float32"), _dict={}, body=None, parent=None):
        """
        
        """
        super(ContactGeometry2DClosedProfile, self).__init__(filename, profile, profile_type="closed", color=color, _dict=_dict, body=body, parent=parent)

    def set_vtk_data(self):
        """

        :return:
        """
        #   points
        self.points = vtk.vtkPoints()
        self.points.SetNumberOfPoints(len(self.profile.nodes))

        #   add points
        for i, node in enumerate(self.profile.nodes):
            if len(node) == 2:
                node = np.append(node, 0.)
            self.points.SetPoint(i, node)

        #   vtkCellArray is a supporting object that explicitly represents cell connectivity.
        #   The cell array structure is a raw integer list of the form:
        #   (n,id1,id2,...,idn, n,id1,id2,...,idn, ...) where n is the number of points in
        #   the cell, and id is a zero-offset index into an associated point list.
        lines = vtk.vtkCellArray()

        lines.InsertNextCell(len(self.profile.nodes) + 1)
        for i in xrange(0, len(self.profile.nodes)):
            lines.InsertCellPoint(i)
        lines.InsertCellPoint(0)

        # vtkPolyData is a data object that is a concrete implementation of vtkDataSet.
        # vtkPolyData represents a geometric structure consisting of vertices, lines,
        # polygons, and/or triangle strips
        self.polygon = vtk.vtkPolyData()
        self.polygon.SetPoints(self.points)
        self.polygon.SetLines(lines)

        polygonMapper = vtk.vtkPolyDataMapper()
        polygonMapper.SetInputData(self.polygon)
        polygonMapper.Update()

        self.vtk_actor = vtk.vtkActor()
        self.vtk_actor.SetMapper(polygonMapper)

        #   scale data
        # self.vtk_actor.SetScale(self.profile.scale)

        #   set color
        self.vtk_actor.GetProperty().SetColor(self.color)

        #   set coordinates
        self.vtk_actor.SetPosition(self.body.R)
        self.vtk_actor.SetOrientation(np.rad2deg(self.body.theta))

class ContactGeometry2DOpenedProfile(ContactGeometry2DProfile):
    """
    classdocs
    """
    def __init__(self, filename, profile, color=np.array([0, 1, 0], dtype="float32"), _dict={}, body=None, parent=None):
        """
        
        """
        super(ContactGeometry2DOpenedProfile, self).__init__(filename, profile, profile_type="opened", color=color, _dict=_dict, body=body, parent=parent)

                 
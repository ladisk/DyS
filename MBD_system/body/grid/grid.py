# coding=utf-8
"""

created by: lskrinjar
date of creation: 16/03/2016
time of creation: 16:34
"""
import itertools
import numpy as np
import vtk
from matplotlib import pyplot as plt


class Grid(object):
    """
    classdocs
    """

    def __init__(self, parent=None):
        """

        :param parent:
        :return:
        """

        self.x_min = 0.0
        self.x_max = +0.6
        self.y_min = -0.6
        self.y_max = +0.6
        self.N = 11

        self.color = np.array([0.8, 0.8, 0.8])

        #   predefined attributes of object
        self.vtk_grid = None
        self.vtk_plane = None
        self.vtk_mapper = None
        self.vtk_actor = None

    def set_vtk_data(self):
        """

        :return:
        """
        div1 = np.linspace(self.x_min, self.x_max, self.N)
        div2 = np.linspace(self.y_min, self.y_max, self.N)
        self.n3 = np.zeros_like(div1)
        self.n1, self.n2 = np.meshgrid(div1, div2)
        self.grid_nodes = np.array([self.n1, self.n2, ])

        #   predefine vtk array
        x = vtk.vtkFloatArray()
        y = vtk.vtkFloatArray()
        z = vtk.vtkFloatArray()

        N = len(self.n1.flatten())

        for n1_i in self.n1.flatten():
            x.InsertNextValue(n1_i)

        for n2_i in self.n2.flatten():
            y.InsertNextValue(n2_i)

        for n3_i in self.n3.flatten():
            z.InsertNextValue(n3_i)

        self.vtk_grid = vtk.vtkRectilinearGrid()
        self.vtk_grid.SetDimensions(N, N, N)
        self.vtk_grid.SetXCoordinates(x)
        self.vtk_grid.SetYCoordinates(y)
        self.vtk_grid.SetZCoordinates(z)

        #   extract a plane from the grid to see what we've got
        self.vtk_plane = vtk.vtkRectilinearGridGeometryFilter()
        self.vtk_plane.SetInputData(self.vtk_grid)
        self.vtk_plane.SetExtent(0, N-1, 0, N-1, 0, 0)

        self.vtk_mapper = vtk.vtkPolyDataMapper()
        self.vtk_mapper.SetInputConnection(self.vtk_plane.GetOutputPort())

        self.vtk_actor = vtk.vtkActor()
        self.vtk_actor.SetMapper(self.vtk_mapper)
        self.vtk_actor.GetProperty().SetRepresentationToWireframe()
        # self.vtk_actor.GetProperty().SetColor(self.color)
        self.vtk_actor.GetProperty().SetColor(.5, .5, .5)

    def update_vtk_data(self):
        """

        :return:
        """

    def plot(self):
        """

        :return:
        """
        fig = plt.figure(figsize=(6, 5),
                     dpi=100,
                     facecolor='w',
                     edgecolor='k')

        plt.plot(grid.n1.flat, grid.n2.flat, ".")
        plt.show()

if __name__ == "__main__":
    grid = Grid()
    # print grid.grid_nodes.shape
    # print grid.grid_nodes


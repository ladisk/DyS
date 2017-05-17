"""

created by: lskrinjar
date of creation: 04/11/2016
time of creation: 11:06
"""
import numpy as np
import vtk
from matplotlib import pyplot as plt


from MBD_system import transform_cs_type
from MBD_system.q_i2R_i_theta_i import q_i2R_i_theta_i
from MBD_system.u_P_lcs2gcs import u_P_lcs2gcs
from MBD_system.Ai_ui_P import Ai_ui_P_vector


class Line(object):
    """
    classdocs
    """
    def __init__(self, node_i=np.array([0, 0], dtype=float), node_j=np.array([1, 0], dtype=float), length=None, direction=np.array([1, 0], dtype=float), color=np.array([0, 1, 0], dtype="float32"), parent=None):
        """
        :param node_i:                node coordinates in body LCS
        :param node_j:                node coordinates in body LCS
        """
        self._parent = parent

        self.node_i = node_i
        self.node_j = node_j
        self.node_list = []
        self.node_list.append(self.node_i)
        self.node_list.append(self.node_j)

        #   length
        self.length = length

        #   direction
        self.direction = direction

        #   length
        # self.length = np.linalg.norm(self.tangent, ord=2)

        #   unit tangent
        # self.tangent_unit = self.tangent / self.length

        #   visualization properties
        self.vtk_actor = None
        self.vtk_mapper = None
        self.color = color
        self.line = None
        self.points = None
        self.line_grid = None

        #   get color from parent
        if hasattr(self._parent, "color"):
            self.color = self._parent.color

        #   line geometry
        self.geometry_nodes = [np.zeros(2), np.zeros(2)]

    def add_attributes_from_dict(self, dictionary):
        """

        :param dictionary:
        :return:
        """
        for key, val in dictionary.iteritems():
            setattr(self, key, val)

        self.update_nodes()

    def update_nodes(self):
        """
        Based on length, the end points of line are reevaluated
        :return:
        """
        angles = [0., np.pi]

        for i, (node, angle) in enumerate(zip(self.node_list, angles)):
            self.node_list[i] = transform_cs_type.transform_polar2cartesian(r=0.5*self.length, phi=angle)#+ self._parent.theta[2]

    def set_vtk_data(self):
        """

        :return:
        """
        #   points
        self.points = vtk.vtkPoints()
        self.points.SetNumberOfPoints(2)
        for i, node in enumerate(self.node_list):
            if len(node) == 2:
                node = np.append(node, 0.)
            self.points.InsertPoint(i, node)

        #   set vtk line object
        self.line = vtk.vtkLine()
        self.line.GetPointIds().SetId(0, 0)
        self.line.GetPointIds().SetId(1, 1)

        self.line_grid = vtk.vtkUnstructuredGrid()
        self.line_grid.Allocate(1, 1)
        self.line_grid.InsertNextCell(self.line.GetCellType(), self.line.GetPointIds())
        self.line_grid.SetPoints(self.points)

        self.vtk_mapper = vtk.vtkDataSetMapper()
        self.vtk_mapper.SetInputData(self.line_grid)

        self.vtk_actor = vtk.vtkActor()
        self.vtk_actor.SetMapper(self.vtk_mapper)
        self.vtk_actor.AddPosition(self._parent.R)
        self.vtk_actor.SetOrientation(np.rad2deg(self._parent.theta))
        self.vtk_actor.GetProperty().SetColor(self.color)

    def update_vtk_data(self):
        """

        :return:
        """

    def evaluate_geometry_nodes(self, q_i):
        """

        :param q_i:     vector of generalized coordinates of a geometry (body)
        :return:
        """
        R_i, theta_i = q_i2R_i_theta_i(q_i)

        for i, uP_lcs in enumerate(self.node_list):
            self.geometry_nodes[i] = R_i + Ai_ui_P_vector(uP_lcs, theta_i)

    def plot_geometry(self, t=None, q_i=None, ax=None, save=False, show=False, color=np.array([0, 0, 0])):
        """

        :param t:
        :param q_i:
        :param ax:
        :param save:
        :param show:
        :param color:
        :return:
        """

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

        #   evaluate geometry
        self.evaluate_geometry_nodes(q_i)

        #   plot geometry
        plt.plot(np.array(self.geometry_nodes)[:, 0], np.array(self.geometry_nodes)[:, 1], color=color)

        if save:
            filename = "t_" + str(t).zfill(2) + ".png"
            fig.savefig( filename)
            fig.clf()


if __name__ == "__main__":
    line = Line()
# coding=utf-8
"""

created by: lskrinjar
date of creation: 25/07/2016
time of creation: 16:12
"""
import itertools
import os


import numpy as np
import vtk


from MBD_system.MBD_system_items import MarkerItem
from MBD_system.q2e_x_jk import q2e_x_jk


class Marker(MarkerItem, vtk.vtkAxesActor):
    """
    classdocs
    """

    _id = itertools.count(0)

    def __init__(self, rP, uP=None, theta0=np.zeros(3), body_id=None, node_id=None, body=None, MBD_item=None, parent=None):
        """

        :return:
        """
        self.marker_id = self._id.next()
        name = "Marker" + str(self.marker_id)
        super(Marker, self).__init__(name, parent=parent)

        #   parent
        self._parent = parent

        #   id
        # self.marker_id = self._id.next()

        #   node id
        self.node_id = node_id

        #   body id
        self.body_id = body_id

        #   position of marker in GCS
        self.rP = rP

        #   position of marker in body LCS
        self.uP = uP

        self.SetConeRadius(0.0)
        self.AxisLabelsOff()

        self._size = 2E-3
        # self.SetShaftTypeToLine()
        self.SetTotalLength(self._size, self._size, self._size)
        # self.SetNormalizedShaftLength(self._size, self._size, 1)
        # self.SetNormalizedTipLength(0., 0., 0.)

        #   pointers
        self.body = body
        self.MBD_item = MBD_item

        #   translate
        transform = vtk.vtkTransform()
        if len(self.rP) == 2:
            self.rP = np.append(self.rP, 0.)

        self.theta0 = theta0
        self.theta = theta0
        # self.theta = self._evaluate_theta()

        transform.Translate(self.rP)
        transform.RotateZ(self.theta[2])
        transform.RotateY(self.theta[1])
        transform.RotateX(self.theta[0])
        self.SetUserTransform(transform)

    def update_vtk_data(self, q, rP=np.zeros(3), theta=np.zeros(3)):
        """

        :return:
        """
        transform = vtk.vtkTransform()

        if len(rP) == 2:
            self.rP[0:2] = rP

        if (theta == np.zeros(3)).all():
            theta = np.rad2deg(self._evaluate_theta(q))
        else:
            theta = theta

        transform.Translate(self.rP)
        transform.RotateZ(theta[2])
        transform.RotateY(theta[1])
        transform.RotateX(theta[0])

        self.SetUserTransform(transform)

    def _evaluate_theta(self, q):
        """

        :return:
        """
        theta = np.zeros(3)
        if self.body is not None:
            if self.body.body_type == "rigid body":
                theta[2] = self.body.theta[2] + self.theta0[2]

            if self.body.body_type == "flexible body":
                e_j_grad = q2e_x_jk(q, self.body.body_id, self.node_id)
                theta[2] = np.arctan2(e_j_grad[1], e_j_grad[0])

        return theta

    def highlightSelected(self):
        """

        :return:
        """
        self.GetXAxisShaftProperty().SetColor(1, 0, 0)
        self.GetYAxisShaftProperty().SetColor(1, 0, 0)
        self.GetZAxisShaftProperty().SetColor(1, 0, 0)

    def unHighlightSelected(self):
        """

        :return:
        """
        self.GetXAxisShaftProperty().SetColor(1, 0, 0)
        self.GetYAxisShaftProperty().SetColor(0, 1, 0)
        self.GetZAxisShaftProperty().SetColor(0, 0, 1)


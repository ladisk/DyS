# coding=utf-8
"""

created by: lskrinjar
date of creation: 20/07/2016
time of creation: 11:36
"""

import datetime
import os
import time
from pprint import pprint
import numpy as np
import vtk
from vtk.qt4.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor
from PyQt4 import QtGui, QtCore, Qt
import sys
from vtk_coordinate_system.vtk_coordinate_system import VTKCoordinateSystem


class VTKWidget(QVTKRenderWindowInteractor):
    """
    classdocs
    Angle units in vtk are always deg
    """

    def __init__(self, MBD_system=None, parent=None):
        """
        Constructor
        """
        super(VTKWidget, self).__init__(parent)

        self._parent = parent
        self._name = "VTKWidget"

        self._typeInfo = "VTK widget"

        #   context menu
        self._parent.setContextMenuPolicy(QtCore.Qt.CustomContextMenu)
        self._parent.customContextMenuRequested.connect(self.contextMenu)

        #   view settings (orthogonal or perspective)
        self.projection = "orthogonal"

        #   info text
        self._show_info = True
        self.text_color = None
        self.textActor_filename = None
        self.textActor_time = None
        self.textActor_step = None
        self.textActor_stepSize = None
        self.textActor_time_date = None
        self.textActor_list = [self.textActor_filename,
                               self.textActor_time,
                               self.textActor_step,
                               self.textActor_stepSize,
                               self.textActor_time_date]

        #   font settings
        self.fontSize = 20

        #   text actors
        self.textactor_filename = None

        # self.vtkWidget = QVTKRenderWindowInteractor(self._parent)

        # self.AddObserver("RightButtonPressEvent", self.contextMenu)

        self.renderer = vtk.vtkRenderer()
        # self.renderer.SetViewport(0.0, 0.0, 0.1, 0.1)
        # print "testing =", self.renderer.GetViewport()
        # print "testing =", self.renderer.GetViewport().GetSize()

        self.GetRenderWindow().AddRenderer(self.renderer)

        self.interactor = self.GetRenderWindow().GetInteractor()
        self.interactor.SetInteractorStyle(vtk.vtkInteractorStyleTrackballCamera())

        self.MBD_system = MBD_system

        #   add ground body for visualization
        self.MBD_system.ground.set_vtk_data()
        if self.MBD_system.ground.grid.vtk_actor is not None:
            self.renderer.AddActor(self.MBD_system.ground.grid.vtk_actor)

        #   orientation GCS
        self.GCS_axesActor = VTKCoordinateSystem()
        self.renderer.AddActor(self.GCS_axesActor)
        
        #   orientation widget
        self.axesActor = vtk.vtkAxesActor()
        self.axesActor.SetScale(4, 4, 4)
        for axis in [self.axesActor.GetXAxisCaptionActor2D(), self.axesActor.GetYAxisCaptionActor2D(), self.axesActor.GetZAxisCaptionActor2D()]:
            axis.GetTextActor().GetTextProperty().ItalicOff()
            axis.GetTextActor().GetTextProperty().SetFontFamilyToCourier()

        self.GCS_orientation = vtk.vtkOrientationMarkerWidget()
        self.GCS_orientation.SetOrientationMarker(self.axesActor)
        self.GCS_orientation.SetInteractor(self.interactor)
        self.GCS_orientation.EnabledOn()
        self.GCS_orientation.InteractiveOff()
        self.GCS_orientation.SetViewport(0, 0, .25, .25)
        self.GCS_orientation.SetInteractor(self.interactor)

        #   set vtk data
        self.set_vtk_data()

        #   camera object
        self.camera = self.renderer.GetActiveCamera()
        if self.projection == 'perspective':
            self.camera.ParallelProjectionOff()
        else:
            self.camera.ParallelProjectionOn()

        self.camera_position = self.camera.GetPosition()
        self.camera_focal_point = self.camera.GetFocalPoint()

        #   set zoom or initial view size
        self._size = 2E-2
        self.camera.SetParallelScale(self._size)
        self.renderer.SetActiveCamera(self.camera)

        colorFunction = vtk.vtkColorTransferFunction()
        colorFunction.SetColorSpaceToHSV()
        colorFunction.HSVWrapOff()
        colorFunction.AddRGBPoint(10, 0.0, 0.0, 1.0)
        colorFunction.AddRGBPoint(20, 1.0, 0.0, 0.0)

        self.scalarBar = vtk.vtkScalarBarActor()
        # self.scalarBar.SeSetScaledText(vtk.vtkTextActor.TEXT_SCALE_MODE_NONE)
        self.scalarBar.SetLookupTable(colorFunction)
        self.scalarBar.SetOrientationToVertical()
        self.scalarBar._dx = 1.
        self.scalarBar._dy = 1.
        self.scalarBar._dx2 = 1.
        self.scalarBar._dy2 = 1.
        self.scalarBar.SetPosition(.9, .1)
        self.scalarBar.SetPosition2(.9, .2)
        self.scalarBar.SetHeight(.8)
        self.scalarBar.SetWidth(0.10)
        self.scalarBar.VisibilityOff()

        #   text property
        self.textProperty = vtk.vtkTextProperty()
        self.textProperty.SetFontFamilyToCourier()
        self.textProperty.SetFontSize(self.fontSize)

        # propT.SetFontFamilyToArial()
        # textProperty.ItalicOff()
        # propT.BoldOn()
        # textProperty.BoldOff()
        self.scalarBar.SetTitleTextProperty(self.textProperty)
        self.scalarBar.SetLabelTextProperty(self.textProperty)
        self.scalarBar.SetLabelFormat("%5.2f")

        # self.scalar_bar_widget = vtk.vtkScalarBarWidget()
        # self.scalar_bar_widget.SetInteractor(self.interactor)
        # self.scalar_bar_widget.SetScalarBarActor(self.scalarBar)
        # self.scalar_bar_widget.SetResizable(False)
        # self.scalar_bar_widget.On()

        self.renderer.AddActor2D(self.scalarBar)

        self.renderer.SetBackground(0, 0, 0)

        self.interactor.Initialize()

        self.displayText()

    def set_vtk_data(self):
        """
        Add actors to renderer
        :return:
        """
        for i, body in enumerate(self.MBD_system.bodies):
            body.set_vtk_data()

            for geom in body.geometry_list:
                if hasattr(geom, "set_vtk_data"):
                    geom.set_vtk_data()

                    if geom.vtk_actor is not None:
                        self.renderer.AddActor(geom.vtk_actor)

            if body.vtk_actor is not None:
                self.renderer.AddActor(body.vtk_actor)

                if body.LCS_actor is not None:
                    self.renderer.AddActor(body.LCS_actor)

            if body.geometry_CS_actor is not None:
                self.renderer.AddActor(body.geometry_CS_actor)

            if body.mesh is not None:
                body.mesh.set_vtk_data()

                if body.mesh.vtk_actor is not None:
                    self.renderer.AddActor(body.mesh.vtk_actor)

                for element in body.mesh.elements:
                    element.set_vtk_data()
                    self.renderer.AddActor(element.vtk_actor)

                    # body.vtkCreateBoxWidget(self.interactor)

        #   create main AABB frame vtk data for visualization
        for contact in self.MBD_system.contacts:
            contact.set_vtk_data(interactor=self.interactor)
            for AABB in contact.AABB_list:
                if AABB is not None:
                    self.renderer.AddActor(AABB.vtk_actor)

            for marker in contact.markers:
                self.renderer.AddActor(marker)

        #   display forces
        for i, force in enumerate(self.MBD_system.forces):
            force.set_vtk_data()

            if force.vtk_actor is not None:
                self.renderer.AddActor(force.vtk_actor)

        #   add spring actors to renderer
        for i, spring in enumerate(self.MBD_system.springs):
            spring.set_vtk_data()

            if spring.vtk_actor is not None:
                self.renderer.AddActor(spring.vtk_actor)

            for marker in spring.markers:
                self.renderer.AddActor(marker)

        #   joints
        for i, joint in enumerate(self.MBD_system.joints):
            joint.set_vtk_data()

            if joint.vtk_actor is not None:
                self.renderer.AddActor(joint.vtk_actor)

            for marker in joint.markers:
                self.renderer.AddActor(marker)

    def update_vtk_data(self, q):
        """
        TODO
        :return:
        """
        print "update_vtk_data()@",__name__

    def contextMenu(self, pos):
        """

        :param pos:
        :return:
        """
        menu = QtGui.QMenu(parent=self)

        refreshAction = menu.addAction("Refresh")
        refreshAction.triggered.connect(self.refresh)
        refreshAction.triggered.connect(lambda: self.MBD_system.update_vtk_data(self.MBD_system.time, self.MBD_system.evaluate_q()))

        restore_q0_Action = menu.addAction("Restore q0")
        restore_q0_Action.triggered.connect(self.MBD_system._restore_initial_conditions)

        menu.addSeparator()

        menu.addMenu(self.viewSubMenu())

        menu.addSeparator()
        if self.GCS_axesActor.GetVisibility():
            _show_GCSAction = QtGui.QAction("Hide GCS", self)
            _show_GCSAction.triggered.connect(self.GCS_axesActor.VisibilityOff)

        else:
            _show_GCSAction = QtGui.QAction("Show GCS", self)
            _show_GCSAction.triggered.connect(self.GCS_axesActor.VisibilityOn)
        menu.addAction(_show_GCSAction)

        menu.addSeparator()
        if self.scalarBar.GetVisibility():
            _show_scalar_bar_Action = QtGui.QAction("Hide scalar bar", self)
            _show_scalar_bar_Action.triggered.connect(self.scalarBar.VisibilityOff)

        else:
            _show_scalar_bar_Action = QtGui.QAction("Show scalar bar", self)
            _show_scalar_bar_Action.triggered.connect(self.scalarBar.VisibilityOn)
        menu.addAction(_show_scalar_bar_Action)

        menu.addSeparator()
        if self.MBD_system.ground.grid.vtk_actor.GetVisibility():
            _show_grid_Action = QtGui.QAction("Hide grid", self)
            _show_grid_Action.triggered.connect(self.MBD_system.ground.grid.vtk_actor.VisibilityOff)

        else:
            _show_grid_Action = QtGui.QAction("Show grid", self)
            _show_grid_Action.triggered.connect(self.MBD_system.ground.grid.vtk_actor.VisibilityOn)
        menu.addAction(_show_grid_Action)

        menu.exec_(self._parent.mapToGlobal(pos))
        self.refresh()

    def viewSubMenu(self):
        """

        :return:
        """
        viewDirectionSubmenu = QtGui.QMenu(self)

        viewDirectionSubmenu.setTitle("View")

        _viewFrontAction = viewDirectionSubmenu.addAction("Front")
        _viewFrontAction.triggered.connect(self._viewFront)

        _viewBackAction = viewDirectionSubmenu.addAction("Back")
        _viewBackAction.triggered.connect(self._viewBack)

        _viewBottomAction = viewDirectionSubmenu.addAction("Bottom")
        _viewBottomAction.triggered.connect(self._viewBottom)

        _viewTopAction = viewDirectionSubmenu.addAction("Top")
        _viewTopAction.triggered.connect(self._viewTop)

        _viewLeftAction = viewDirectionSubmenu.addAction("Left")
        _viewLeftAction.triggered.connect(self._viewLeft)

        _viewRightAction = viewDirectionSubmenu.addAction("Right")
        _viewRightAction.triggered.connect(self._viewRight)

        viewDirectionSubmenu.addSeparator()
        _viewIsometricAction = viewDirectionSubmenu.addAction("Isometric")
        _viewIsometricAction.triggered.connect(self._viewIsometric)

        viewDirectionSubmenu.addSeparator()
        _refitAction = viewDirectionSubmenu.addAction("Refit")
        _refitAction.triggered.connect(self._refit)

        return viewDirectionSubmenu

    def _viewFront(self):
        """

        :return:
        """
        self.camera.SetPosition(self.camera_position)
        self.camera.SetFocalPoint(self.camera_focal_point)
        self.camera.SetViewUp(np.array([0, 1, 0], dtype="float"))
        self.camera.ParallelProjectionOn()
        self.camera.SetParallelScale(self._size)

    def _viewBack(self):
        """

        :return:
        """
        self.camera.SetPosition(np.array([0, 0, -1], dtype="float"))
        self.camera.SetFocalPoint(self.camera_focal_point)
        self.camera.SetViewUp(np.array([0, 1, 0], dtype="float"))
        self.camera.ParallelProjectionOn()
        self.camera.SetParallelScale(self._size)

    def _viewBottom(self):
        """

        :return:
        """
        self.camera.SetPosition(np.array([0, -1, 0], dtype="float"))
        self.camera.SetFocalPoint(self.camera_focal_point)
        self.camera.SetViewUp(np.array([0, 0, 1], dtype="float"))
        self.camera.ParallelProjectionOn()
        self.camera.SetParallelScale(self._size)

    def _viewTop(self):
        """

        :return:
        """
        self.camera.SetPosition(np.array([0, +1, 0], dtype="float"))
        self.camera.SetFocalPoint(self.camera_focal_point)
        self.camera.SetViewUp(np.array([0, 0, -1], dtype="float"))
        self.camera.ParallelProjectionOn()
        self.camera.SetParallelScale(self._size)

    def _viewLeft(self):
        """

        :return:
        """
        self.camera.SetPosition(np.array([-1, 0, 0], dtype="float"))
        self.camera.SetFocalPoint(self.camera_focal_point)
        self.camera.SetViewUp(np.array([0, +1, 0], dtype="float"))
        self.camera.ParallelProjectionOn()
        self.camera.SetParallelScale(self._size)

    def _viewRight(self):
        """

        :return:
        """
        self.camera.SetPosition(np.array([+1, 0, 0], dtype="float"))
        self.camera.SetFocalPoint(self.camera_focal_point)
        self.camera.SetViewUp(np.array([0, +1, 0], dtype="float"))
        self.camera.ParallelProjectionOn()
        self.camera.SetParallelScale(self._size)

    def _viewIsometric(self):
        """

        :return:
        """
        self.camera.SetPosition(np.array([+1, +1, +1], dtype="float"))
        self.camera.SetFocalPoint(self.camera_focal_point)
        self.camera.SetViewUp(np.array([0, +1, 0], dtype="float"))
        self.camera.ParallelProjectionOn()
        self.camera.SetParallelScale(self._size)

    def _saveSnapShot(self):
        """

        :return:
        """
        # self.renderer.SetBackground(1, 1, 1)

        windowToImageFilter = vtk.vtkWindowToImageFilter()
        windowToImageFilter.SetInput(self.GetRenderWindow())
        # windowToImageFilter.SetMagnification(3)
        windowToImageFilter.SetInputBufferTypeToRGBA()# also record the alpha (transparency) channel
        windowToImageFilter.ReadFrontBufferOff()
        windowToImageFilter.Update()

        filename = "screenshot_" + datetime.datetime.now().strftime('%Y.%m.%d_%H.%M.%S') + ".jpg"

        writer = vtk.vtkJPEGWriter()
        writer.SetFileName(filename)
        writer.SetInputConnection(windowToImageFilter.GetOutputPort())
        writer.Write()

        print "Snapshot saved to file %s"%os.path.abspath(filename)

    def _refit(self):
        """

        :return:
        """
        print "TODO!"

    def resizeEvent(self, QResizeEvent):
        """

        :param QResizeEvent:
        :return:
        """
        if self.textActor_filename is not None:
            self.textActor_filename.SetPosition(self.textActor_filename.dx0, self.geometry().height() - self.textActor_filename.dy0)

        if self.textActor_time is not None:
            self.textActor_time.SetPosition(self.textActor_time.dx0, self.geometry().height() - self.textActor_time.dy0)

        if self.textActor_step is not None:
            self.textActor_step.SetPosition(self.textActor_step.dx0, self.geometry().height() - self.textActor_step.dy0)

        if self.textActor_time_date is not None:
            self.textActor_time_date.SetPosition(self.geometry().width() - self.textActor_time_date.dx0, self.geometry().height() - self.textActor_time_date.dy0)

        if self.textActor_stepSize is not None:
            self.textActor_stepSize.SetPosition(self.textActor_stepSize.dx0, self.geometry().height() - self.textActor_stepSize.dy0)

        self.textProperty.SetFontSize(self.fontSize)
        # self.scalarBar.SetTitleTextProperty(self.textProperty)
        # self.scalarBar.SetLabelTextProperty(self.textProperty)

        self.scalarBar.SetPosition(.9, .1)
        self.scalarBar.SetPosition2(.9, .2)
        self.scalarBar.SetHeight(.8)
        self.scalarBar.SetWidth(0.10)

    def displayText(self):
        """

        :return:
        """
        dx0 = 10.

        #   filename
        self.textActor_filename = vtk.vtkTextActor()
        self.textActor_filename.SetInput(self.MBD_system._name)

        textProperty = self.textActor_filename.GetTextProperty()
        textProperty.SetFontFamilyToCourier()
        textProperty.SetFontSize(self.fontSize)

        self.textActor_filename.dx0 = dx0
        self.textActor_filename.dy0 = 30.

        self.textActor_filename.SetPosition(self.textActor_filename.dx0, self.geometry().height() - self.textActor_filename.dy0)
        self.renderer.AddActor(self.textActor_filename)

        #   simulation time
        self.textActor_time = vtk.vtkTextActor()
        self.textActor_time.SetInput("t: " + str(self.MBD_system.time))

        textProperty = self.textActor_time.GetTextProperty()
        textProperty.SetFontFamilyToCourier()
        textProperty.SetFontSize(self.fontSize)

        self.textActor_time.dx0 = dx0
        self.textActor_time.dy0 = 60.
        self.textActor_time.SetPosition(self.textActor_time.dx0, self.geometry().height() - self.textActor_time.dy0)
        self.renderer.AddActor(self.textActor_time)

        #   step
        self.textActor_step = vtk.vtkTextActor()
        self.textActor_step.SetInput("step: " + str(self.MBD_system.step_num))

        textProperty = self.textActor_step.GetTextProperty()
        textProperty.SetFontFamilyToCourier()
        textProperty.SetFontSize(self.fontSize)

        self.textActor_step.dx0 = dx0
        self.textActor_step.dy0 = 90.
        self.textActor_step.SetPosition(self.textActor_step.dx0, self.geometry().height() - self.textActor_step.dy0)
        self.renderer.AddActor(self.textActor_step)

        #   step size
        self.textActor_stepSize = vtk.vtkTextActor()
        self.textActor_stepSize.SetInput("h: " + str(self.MBD_system.h))

        textProperty = self.textActor_stepSize.GetTextProperty()
        textProperty.SetFontFamilyToCourier()
        textProperty.SetFontSize(self.fontSize)

        self.textActor_stepSize.dx0 = dx0
        self.textActor_stepSize.dy0 = 120.
        self.textActor_stepSize.SetPosition(self.textActor_stepSize.dx0, self.geometry().height() - self.textActor_stepSize.dy0)
        self.renderer.AddActor(self.textActor_stepSize)

        #   time and date
        self.textActor_time_date = vtk.vtkTextActor()
        self.textActor_time_date.SetInput(str(time.strftime("%b %d %Y %H:%M:%S")))

        textProperty = self.textActor_time_date.GetTextProperty()
        textProperty.SetFontFamilyToCourier()
        textProperty.SetFontSize(self.fontSize)

        self.textActor_time_date.dx0 = 300.
        self.textActor_time_date.dy0 = 30.
        self.textActor_time_date.SetPosition(self.geometry().width() - self.textActor_time_date.dx0, self.geometry().height() - self.textActor_time_date.dy0)
        self.renderer.AddActor(self.textActor_time_date)

    def refresh(self, step=None, h=None):
        """

        :return:
        """
        # q = self.MBD_system.evaluate_q()
        # self.MBD_system.update_vtk_data(self.MBD_system.time, q)

        self.textActor_time_date.SetInput(str(time.strftime("%b %d %Y %H:%M:%S")))

        #   simulation time
        self.textActor_time.SetInput("t: " + str(self.MBD_system.time))

        #   step
        if step is None:
            step = self.MBD_system.step_num
        self.textActor_step.SetInput("step: " + str(step))

        #   step size
        if h is None:
            h = self.MBD_system.h

        self.textActor_stepSize.SetInput("h: " + str(h))

        self.update()

    def setTextColor(self, color):
        """
        
        :param color: 
        :return: 
        """
        for text_actor in [self.textActor_filename,
                           self.textActor_time,
                           self.textActor_step,
                           self.textActor_stepSize,
                           self.textActor_time_date]:
            text_actor.GetTextProperty().SetColor(color)

        #   change color of coordinate system labels
        for axis in [self.axesActor.GetXAxisCaptionActor2D(), self.axesActor.GetYAxisCaptionActor2D(), self.axesActor.GetZAxisCaptionActor2D()]:
            axis.GetTextActor().GetTextProperty().SetColor(color)

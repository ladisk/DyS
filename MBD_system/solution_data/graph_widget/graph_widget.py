"""
Created on 22. apr. 2016

@author: luka.skrinjar
"""
import sys
from pprint import pprint
import numpy as np
import scipy.ndimage as ndi
import itertools
import vtk
from vtk.qt4.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor
from PyQt4 import QtGui, QtCore


from graph_widget_ui import Ui_Form

try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    def _fromUtf8(s):
        return s


class GraphWidget(QVTKRenderWindowInteractor):#QtGui.QMainWindow, GraphicsWindow, GraphicsView, QVTKRenderWindowInteractor, QtGui.QWidget
    """
    classdocs
    """
    __id = itertools.count(0)

    def __init__(self, measure_obj=None, parent=None):
        """
        Constructor
        """
        super(GraphWidget, self).__init__(parent=parent)
        
        #   parent
        self._parent = parent

        # self.frame = QtGui.QFrame()
        #
        # self.vl = QtGui.QVBoxLayout()
        self.view = vtk.vtkContextView()
        self.view.GetRenderer().SetBackground(1.0, 1.0, 1.0)
        # self.vtkWidget = QVTKRenderWindowInteractor(self.frame)
        #
        # self.vl.addWidget(self.vtkWidget)

        self.view.SetInteractor(self.GetRenderWindow().GetInteractor())
        self.SetRenderWindow(self.view.GetRenderWindow())

        # self.ren = vtk.vtkRenderer()
        # self.vtkWidget.GetRenderWindow().AddRenderer(self.ren)
        # self.interactor = self.vtkWidget.GetRenderWindow().GetInteractor()



        # self.frame.setLayout(self.vl)
        # self.setCentralWidget(self.frame)

        self.xy_data = vtk.vtkTable()
        self.arrX = vtk.vtkFloatArray()
        self.arrX.SetName("X Axis")
        self.xy_data.AddColumn(self.arrX)

        self.arrY = vtk.vtkFloatArray()
        self.arrY.SetName("Y value")
        self.xy_data.AddColumn(self.arrY)

        numPoints = 20
        inc = 7.5 / (numPoints - 1)
        self.xy_data.SetNumberOfRows(numPoints)
        for i in range(numPoints):
            self.xy_data.SetValue(i, 0, i * inc)
            self.xy_data.SetValue(i, 1, i * inc)

        # self.view = vtk.vtkContextView()
        # self.view.GetRenderer().SetBackground(1.0, 0.0, 1.0)
        #
        self.chart = vtk.vtkChartXY()
        self.chart.GetAxis(0).SetGridVisible(False)
        self.chart.GetAxis(1).SetGridVisible(False)

        self.view.GetScene().AddItem(self.chart)
        self.line = self.chart.AddPlot(vtk.vtkChart.LINE)
        self.line.SetInputData(self.xy_data, 0, 1)
        #
        # self.ren = self.view.GetRenderer()
        # self.renWin = self.view.GetRenderWindow()
        # self.renWin.AddRenderer(self.ren)

        # iren = self.view.GetInteractor()
        # iren.SetRenderWindow(self.view.GetRenderWindow())
        #
        # iren.Initialize()
        # self.view.GetRenderWindow().Render()


if __name__ == "__main__":
    app = QtGui.QApplication(sys.argv)
    win = GraphWidget()
    win.x = np.linspace(0, 100, 100)
    win.y = np.random.random(size=100)
    # win.ui.graphicsView.plot(win.x, win.y)
    # win.show()
    win.view.GetInteractor().Start()
    # win.interactor.Initialize()
    sys.exit(app.exec_())

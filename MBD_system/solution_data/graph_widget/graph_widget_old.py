"""
Created on 22. apr. 2016

@author: luka.skrinjar
"""
import sys
from pprint import pprint
import numpy as np
import scipy.ndimage as ndi
import itertools
import pyqtgraph as pg
from PyQt4 import QtGui, QtCore


from graph_widget_ui import Ui_Form


class GraphWidget(QtGui.QWidget):#GraphicsWindow, GraphicsView
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
        
        #    id
        self._graph_id = self.__id.next()

        #   backgrund color
        pg.setConfigOption('background', 'w')

        #    init ui file
        self.ui = Ui_Form()
        self.ui.setupUi(self)
        
        #   set line color, width
        self.pen = pg.mkPen(width=1.5, color=(200, 200, 255))
        self.pen_crosshair = pg.mkPen(width=1., color=(200, 200, 200))
        
        #    crosshair
        self.vLine = pg.InfiniteLine(angle=90, pen=self.pen_crosshair, movable=False)
        self.hLine = pg.InfiniteLine(angle=0, pen=self.pen_crosshair, movable=False)
        self.ui.graphicsView.addItem(self.vLine, ignoreBounds=True)
        self.ui.graphicsView.addItem(self.hLine, ignoreBounds=True)
        
        #   create label
        self.label_anchor = (0, 1)
#         self.label = pg.LabelItem(justify='right')
        
        self.label = pg.TextItem("x= , y= ", anchor=self.label_anchor)
        self.ui.graphicsView.addItem(self.label)
        
        #   set size
        _screen_geometry = QtGui.QDesktopWidget().screenGeometry()
        self.resize(int(.2*_screen_geometry.width()), int(.2*_screen_geometry.height()))
 
        #   disable context menu
        self.ui.graphicsView.setMenuEnabled(enableMenu=False)
        #   disable mouse for axes x and y
        self.ui.graphicsView.setMouseEnabled(False, False)

        # self.ui.graphicsView.setConfigOption('Mmiddlebuttonzoom', False)
 
        # self.graph_widget.ui.graphicsView._menuEnabled = False
        self.ui.graphicsView.hideButtons()

        #   show all axis
        self.ui.graphicsView.showGrid(x=True, y=True, alpha=.5)
        #   axis label
        self.x_axis = self.ui.graphicsView.getAxis("bottom")
        self.labelStyle = {'font-size': '12pt'}
        #   axis initial range
        # self.ui.graphicsView.setYRange(-1, +1, padding=0)

        #   move
        self.move(self._graph_id * self.width(), _screen_geometry.height() - self.height() - 100)

        #   visibility of widget
        self._visible = True

        #   x, y coordinates at mouse position
        self.ui.graphicsView.scene().sigMouseMoved.connect(self.mouseCoordinates)
        
        #    data
        self.x = []
        self.y = []
    
        #    measure object
        self.measure_obj = measure_obj

    def mouseCoordinates(self, pos):
        """

        :return:
        """
        if self.ui.graphicsView.sceneBoundingRect().contains(pos):
            mousePoint = self.ui.graphicsView.centralWidget.vb.mapSceneToView(pos)

            index = int(mousePoint.x())
            if self.measure_obj is None:
                if index > 0 and index < len(self.x):
                    self.label.setText("x=%0.1f, y=%.6E" % (self.x[index], self.y[index]))
                    
                    self.vLine.setPos(mousePoint.x())
                    self.hLine.setPos(mousePoint.y())
            
            else:
                if index > 0 and index < len(self.measure_obj.x):
                    self.label.setText("x=%0.1f, y=%.6E" % (self.measure_obj.x[index], self.measure_obj.y[index]))

                    self.vLine.setPos(mousePoint.x())
                    self.hLine.setPos(mousePoint.y())

    def close(self, event):
        """

        :param event:
        :return:
        """
        self._visible = False
        
        
if __name__ == "__main__":
    app = QtGui.QApplication(sys.argv)
    win = GraphWidget()
    win.x = np.linspace(0, 100, 100)
    win.y = np.random.random(size=100)
    win.ui.graphicsView.plot(win.x, win.y)
    win.show()
    sys.exit(app.exec_())

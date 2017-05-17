"""

created by: lskrinjar
date of creation: 27/01/2016
time of creation: 12:41
"""
from pprint import pprint
import itertools
from PyQt4 import QtCore, QtGui


from MBD_system.solution_data.graph_widget.graph_widget import GraphWidget
from MBD_system.MBD_system_items import MeasureItem


class Measure(MeasureItem):
    """
    classdocs
    """
    __id = itertools.count(0)

    def __init__(self, y_variable, x_variable="time", body_id_i=None, body_id_j=None, name=None, parent=None):
        """
        Constructor
        """
        super(Measure, self).__init__(name, parent=parent)

        #   parent
        self._parent = parent

        #   id
        self.measure_id = self.__id.next()

        #   name
        self._name = name
        
        #    id of item to measure its property value
        self.body_id_i = body_id_i
        self.body_id_j = body_id_j

        #    variable to measure
        self.x_variable = x_variable
        self.y_variable = y_variable

        #   init plot containers
        #   x axis (default is time)
        self.x = [0]
        #   y axis is defined by user
        self.y = []

        #   graph widget
        # print "self._parent._parent.main_app =", self._parent._parent.main_app
        # print "self._parent._parent.main_app =", self._parent._parent.main_app, type(self._parent._parent.main_app)
        # pprint(vars(self._parent._parent.main_app))
        self.graph_widget = GraphWidget(measure_obj=self, parent=self._parent._parent.dys)#self._parent._parent.dys
        self.graph_widget.setWindowTitle("Measure_" + str(self.measure_id))

        #   change axis label
        self.graph_widget.x_axis.setLabel(text=self.x_variable, units=None, unitPrefix=None)

        #   get y axis as object
        self.y_axis = self.graph_widget.ui.graphicsView.getAxis("left")

    def _track_data(self, step, h, t, q):
        """

        :param step:
        :param h:
        :param t:
        :param q:
        :return:
        """
        #   this method is created for every type of MBD object to measure its variables
        self._measure(step, h, t, q)
    
    def _paintGL(self):
        """
        
        """
        if self.graph_widget._visible:
            try:
                self.graph_widget.ui.graphicsView.plot(self.x[0:len(self.y)], self.y, pen=self.graph_widget.pen, clear=False)
            except:
#                 pass
#                 print "vectors X and Y not equal size, ", self._name
#                 print "x =", len(self.x)
#                 print "y =", len(self.y)
                raise Warning, "vectors X and Y not equal size, %s" %self._name
    #
    def reset(self):
        """
        Function resets specific attributes related to simulation
        :return:
        """
        self.x = []
        self.y = []

        self.graph_widget.ui.graphicsView.clear()

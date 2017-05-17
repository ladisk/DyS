"""
Created on 28. apr. 2016

@author: luka.skrinjar
"""
import time
import numpy as np


from measure import Measure


class MeasureContact(Measure):
    """
    classdocs
    """
    
    def __init__(self, y_variable, contact_id, x_variable="time", name=None, parent=None):
        """
        Constructor
        """
        super(MeasureContact, self).__init__(y_variable, x_variable=x_variable, parent=parent)
        
        #   set name
        if name is None:
            self._name = "contact_" + str(contact_id) + "_" + y_variable
        else:
            self._name = name
        
        #   contact id
        self.contact_id = contact_id
        
        #   supported types
        self.y_variables = ["Fx",
                            "Fy",
                              "Fn_i",
                              "Ft_i",
                              "Fn_j",
                              "Ft_j",
                              "v_n",
                              "v_t",
                              "delta",
                              "status"]

        #    pointer to contact object
        self.contact = self._parent._parent.contacts[self.contact_id]

        #   set y axis label and units
        y_axis = self.graph_widget.ui.graphicsView.getAxis("left")
        if "F" in self.y_variable:
            y_axis.setLabel(text="F", units="N", unitPrefix=None, **self.graph_widget.labelStyle)

            if self.y_variable[1] in ["x", "y"]:
                y_axis.setLabel(text="F"+self.y_variable[1], units="N", unitPrefix=None, **self.graph_widget.labelStyle)

        if self.y_variable == "v_n":
            y_axis.setLabel(text="v_n", units="m/s", unitPrefix=None, **self.graph_widget.labelStyle)

        if self.y_variable == "v_t":
            y_axis.setLabel(text="v_t", units="m/s", unitPrefix=None, **self.graph_widget.labelStyle)

        if "delta" in self.y_variable:
            y_axis.setLabel(text="delta", units="m", unitPrefix=None, **self.graph_widget.labelStyle)
            
        if "status" in self.y_variable:
            y_axis.setLabel(text="status", units=" ", unitPrefix=None, **self.graph_widget.labelStyle)
        
    def _measure(self, step, h, t, q):
        """

        :param t:
        :return q:  vector of state of MBD system
        """
        #   x data to plot
        if self.x_variable == "time":
            self.x.append(t)
        elif self.x_variable == "step":
            self.x.append(step)

        #   y component to plot
        if self.y_variable == "Fx":
            self.y = np.vstack(self.contact.solution_data._F_solution_container)[:, 0]

        if self.y_variable == "Fy":
            self.y = np.vstack(self.contact.solution_data._F_solution_container)[:, 1]

        if self.y_variable == "Fn_i":
            self.y = np.array(self.contact._Fn_solution_container)

        if self.y_variable == "delta":
            self.y = np.array(self.contact._distance_solution_container)
        
        if self.y_variable == "status":
            self.y = np.array(self.contact._status_container)
        
        if self.y_variable == "v_n":
            self.y = np.array(self.contact._dqn_solution_container)

        if self.y_variable == "v_t":
            self.y = np.array(self.contact._dqt_solution_container)
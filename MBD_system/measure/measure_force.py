"""
Created on 28. apr. 2016

@author: luka.skrinjar
"""

import numpy as np


from measure import Measure


class MeasureForce(Measure):
    """
    classdocs
    """
    
    def __init__(self, y_variable, contact_id, x_variable="time", name=None, parent=None):
        """
        Constructor
        """
        super(MeasureForce, self).__init__(y_variable, parent=parent)
        
        #   set name
        if name is None:
            self._name = "force_" + str(contact_id) + "_" + y_variable
        else:
            self._name = name
        
        #   force id
        self.force_id = contact_id
        
        #   supported types
        self.y_variables = ["F",
                              "Fx",
                              "Fy"]
        
        #    pointer to contact object
        self.force = self._parent._parent.forces[self.force_id]
        
    def _measure(self,step, h, t, q):
        """

        :param t:
        :return q:  vector of state of MBD system
        """
        self.x = self.force._t_solution_container
        #   y component to plot
        if self.y_variable == "Fx":
            self.y = np.array(self.force._Fx_solution_container)

            
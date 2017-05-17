"""
Created on 28. apr. 2016

@author: luka.skrinjar
"""

from measure import Measure
from MBD_system.q2qdq_i import q2qdq_i


class MeasureBody(Measure):
    """
    classdocs
    """
    
    def __init__(self, y_variable, body_id, x_variable="time", node_id=None, name=None, parent=None):
        """
        Constructor
        """
        super(MeasureBody, self).__init__(y_variable, parent=parent)

        #   set name
        if name is None:
            self._name = "body_" + str(body_id) + "_" + y_variable
        else:
            self._name = name

        #   body id
        self.body_id = body_id

        #   node is
        self.node_id = node_id

        #   body pointer
        self.body = self._parent._parent.bodies[body_id]

        #   supported types
        self.y_variables = ["Rx",
                              "Ry",
                              "theta",
                              "dRx",
                              "dRy",
                              "dtheta",
                              "potential_energy",
                              "kinetic_energy",
                              "translational kinetic energy",
                              "angular kinetic energy",
                              "mechanical_energy",
                              "translational momentum",
                              "angular momentum"]

        if self.y_variable == "Rx":
            self._y_index = 0
            self.y_axis.setLabel(text="Rx", units="m", unitPrefix=None, **self.graph_widget.labelStyle)

        if self.y_variable == "Ry":
            self._y_index = 1
            self.y_axis.setLabel(text="Ry", units="m", unitPrefix=None, **self.graph_widget.labelStyle)

        if self.y_variable == "theta":
            self._y_index = 2
            self.y_axis.setLabel(text="theta", units="rad", unitPrefix=None, **self.graph_widget.labelStyle)

        if self.y_variable == "dRx":
            self._y_index = 3
            self.y_axis.setLabel(text="dRx", units="m/s", unitPrefix=None, **self.graph_widget.labelStyle)

        if self.y_variable == "dRy":
            self._y_index = 4
            self.y_axis.setLabel(text="dRy", units="m/s", unitPrefix=None, **self.graph_widget.labelStyle)

        #   plot y axis label
        if self.y_variable == "kinetic_energy":
            self.y_axis.setLabel(text="Ek", units="J", unitPrefix=None, **self.graph_widget.labelStyle)

        #   plot y axis label
        if self.y_variable == "mechanical_energy":
            self.y_axis.setLabel(text="Em", units="J", unitPrefix=None, **self.graph_widget.labelStyle)

    def _measure(self, step, h, t, q):
        """

        :param t:
        :return q:  vector of state of MBD system
        """
        if self.x_variable == "time":
            self.x.append(t)
        elif self.x_variable == "step":
            self.x.append(step)

        #   body q_i
        q_i, dq_i = q2qdq_i(q, self.body_id)

        #   y component to plot
        if self.y_variable == "Rx" or self.y_variable == "Ry" or self.y_variable == "theta" or self.y_variable == "dRx" or self.y_variable == "dRy" or self.y_variable == "dtheta":
            #   rigid body
            if self.node_id is None:
                y_i = q_i[self._y_index]
            #   flexible body - node id from mesh
            else:
                self.body.mesh.e_j(node_id=self.node_id, e=q_i)
                y_i = q_i[self._y_index]

        elif self.y_variable == "kinetic_energy":
            y_i = self.body._kinetic_energy
        
        elif self.y_variable == "mechanical_energy":
            y_i = self.body.evaluate_mechanical_energy()

        else:
            y_i = None
            print "error in measure_body.py"

        #   add to list (container)
        self.y.append(y_i)

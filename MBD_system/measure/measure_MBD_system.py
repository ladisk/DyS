"""

created by: lskrinjar
date of creation: 29/04/2016
time of creation: 19:23
"""

from measure import Measure


class MeasureMBDsystem(Measure):
    """
    classdocs
    """

    def __init__(self, measure_type, MBD_system=None, name=None, parent=None):
        """
        Constructor
        """
        super(MeasureMBDsystem, self).__init__(measure_type, parent=parent)

        #   set name
        if name is None:
            self._name = "MBD_system_" + measure_type
        else:
            self._name = name

        #   pointer to MBD system object
        if MBD_system is None:
            self.MBD_system = self._parent._parent
        else:
            self.MBD_system = MBD_system

        #   over ride initial settings for variable x
        self.x = []

        #   supported types of measurements
        self.y_variables = ["kinetic_energy",
                              "potential_energy",
                              "mechanical_energy"
                              "translational momentum",
                              "angular momentum"]

        y_axis = self.graph_widget.ui.graphicsView.getAxis("left")
        if self.y_variable == "mechanical_energy":
            y_axis.setLabel(text="Em", units="J", unitPrefix=None, **self.graph_widget.labelStyle)

    def _measure(self, step, h, t, q):
        """

        :param t:
        :return q:  vector of state of MBD system
        """
        if self.x_variable == "time":
            self.x.append(t)
        elif self.x_variable == "step":
            self.x.append(step)

        #   y component to plot
        if self.y_variable == "kinetic_energy":
            y_i = self.MBD_system._kinetic_energy

        if self.y_variable == "mechanical_energy":
            y_i = self.MBD_system._mechanical_energy

        #   add to list (container)
        self.y.append(y_i)
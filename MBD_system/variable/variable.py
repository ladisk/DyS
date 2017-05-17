"""
Created on 8. jun. 2016

@author: luka.skrinjar
"""
import numpy as np
import matplotlib.pyplot as plt


from MBD_system.MBD_system_items import VariableItem


class Variable(VariableItem):
    """
    classdocs
    """


    def __init__(self, MBD_item_type, item_id, variable, x_n, x_l, x_u, parent=None):
        """
        Constructor
        """
        super(Variable, self).__init__(variable+"_"+MBD_item_type+"_"+str(item_id), parent)
        #   parent
        self._parent = parent

        #   type of object
        self.MBD_item_type = MBD_item_type

        #   id
        self.item_id = item_id

        #   variable name
        self.variable = variable

        #   values for distribution
        #   nominal value
        self.x_n = x_n

        #   nominal value
        self.x_l = x_l

        #   nominal value
        self.x_u = x_u

        #   standard deviation
        self.standard_deviation = self.x_u - self.x_l

        #   pointer to the MBD item
        self.MBD_item = None
        if self.MBD_item_type == "body":
            self.MBD_item = self._parent._parent.bodies[self.item_id]

        elif self.MBD_item_type == "force":
            self.MBD_item = self._parent._parent.forces[self.item_id]

        elif self.MBD_item_type == "joint":
            self.MBD_item = self._parent._parent.joints[self.item_id]

        elif self.MBD_item_type == "contact":
            self.MBD_item = self._parent._parent.contact[self.item_id]

        elif self.MBD_item_type == "spring":
            self.MBD_item = self._parent._parent.springs[self.item_id]

        elif self.MBD_item_type == "motion":
            self.MBD_item = self._parent._parent.motions[self.item_id]

        else:
            print "MBD_item not defined!"

    def set_value(self):
        """

        :return:
        """

        self.value = np.random.normal(loc=self.x_n, scale=self.standard_deviation, size=None)

        #   assign value to MBD item object
        if self.MBD_item is not None:
            setattr(self.MBD_item, self.variable, self.value)


if __name__ == "__main__":
    n = 100

    var = Variable(MBD_item_type="Undefined",
                   item_id=0,
                   variable="test",
                   x_n=10,
                   x_l=8,
                   x_u=12,
                   parent=None)

    values = []
    for i in range(0, n):
        var.set_value()

        values.append(var.value)

    hist, bin_edges = np.histogram(values, density=True)

    plt.hist(values, bins='auto')  # plt.hist passes it's arguments to np.histogram
    plt.title("Histogram with 'auto' bins")
    plt.show()
__author__ = 'lskrinjar'

from MBD_system.contact.contact import Contact

class PinSlotClearanceJointLinear(Contact):
    """

    """
    def __init__(self, _type, body_id_i, body_id_j, u_iP, u_jP, u_jR, R0_i, h0_jP, , h0_jR=None, parent=None):
        """

        :param _type:
        :param body_id_i:
        :param body_id_j:
        :return:
        """
        super(Contact, self).__init__(_type, body_id_i, body_id_j, name=self._name, properties_dict=properties_dict, parent=parent)

        self._parent = parent

        #   body i (pin) geometry parameters
        self.u_iP = u_iP
        #   radius of a pin
        self.R0_i = R0_i

        #   body j (slot) geometry data
        #   vector to first point of center of a slot
        self.u_jP = u_jP
        #   vector to second point of center of a slot
        self.u_jR = u_jR
        #   width of a slot at first point of center of a slot
        self.h0_jP = h0_jP
        #   width of a slot at second point of center of a slot
        if h0_jR is None:
            self.h0_jR = h0_jP
        else:
            self.h0_jR = h0_jR


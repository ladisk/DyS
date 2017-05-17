"""

created by: lskrinjar
date of creation: 09/08/2016
time of creation: 14:56
"""


class GlobalVariables(object):
    """

    Stores the global variables, accessible from all modules.

    Attr:
        None
    """

    def __init__(self):
        """
        Initiates the GlobalVariables module class.
        """
        #   folder path
        self.MBDsystem_folder = ""

        #   vector of coordinates per MBD body/finite element
        self.q_i_dim = []

        #   vector of constrained coordinates per kinematic constraint (joint, support) [rows, cols]
        self.C_q_i_dim = []

        #   step number of numerical simulation
        self.step = None

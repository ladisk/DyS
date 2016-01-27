"""
Created on 27. jan. 2016

@author: lskrinjar
"""
import logging
import itertools
import time
from operator import attrgetter
from pprint import pprint
import numpy as np
from matplotlib import pyplot as plt

from MBD_system.contact.distance.distance import Distance
from MBD_system.force.force import Force
from MBD_system.q2R_i import q2R_i
from MBD_system.q2theta_i import q2theta_i
from MBD_system.q2dR_i import q2dR_i
from MBD_system.q2dtheta_i import q2dtheta_i
from MBD_system.q2q_body import q2q_body
from MBD_system.transform_cs import cm_lcs2gcs, gcs2cm_lcs
from MBD_system.dr_contact_point_uP import dr_contact_point_uP
from MBD_system.contact_model.contact_model import ContactModel
from MBD_system.contact_model.contact_model_cylinder import ContactModelCylinder
from MBD_system.contact.contact import Contact
from MBD_system.A import A_matrix
from MBD_system.Ai_ui_P import Ai_ui_P_vector
from MBD_system.transform_cs import uP_gcs2lcs
from simulation_control_widget.opengl_widget.marker.marker import Marker
from MBD_system.contact.distance.distance_revolute_clearance_joint import DistanceRevoluteClearanceJoint

class GeneralContact(Contact):
    """
    classdocs
    """
    __id = itertools.count(0)

    def __init__(self, _type, body_id_i, body_id_j, properties_dict=[], parent=None):
        """
        Constructor of class contact of revolute clearance joint
        :param _type:       type of clearance joint or contact
        :param body_id_i:   id of hole body
        :param body_id_j:   id of pin body
        :param properties_dict: additioanl parameters to override default values or add new parameters
        """
        #   parent
        self._parent = parent

        #   contact point in GCS
        self.u_P_GCS = None

        #   AABB - Axis Aligned Bounding Box properties
        #   initialize a list of AABB
        self.AABB_i = None
        self.AABB_j = None
        self.AABB_list = []
        #   initialize a list of contact pairs
        self.AABB_list_of_overlap_pairs = []
        #   contact status attributes
        self.AABB_overlay_detected = False

        #   reset to initial value
        self.reset()

    def reset(self):
        """
        Function resets object attributes to initial value to be used before simulation
        :return:
        """


if __name__ == "__main__":

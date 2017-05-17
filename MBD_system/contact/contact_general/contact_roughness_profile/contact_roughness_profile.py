"""

created by: lskrinjar
date of creation: 28/03/2016
time of creation: 13:11
"""
from pprint import pprint
import numpy as np

from MBD_system.contact.contact_general.contact_general import GeneralContact
from MBD_system.contact.contact_general.contact_roughness_profile.roughness_profile.roughness_profile import RoughnessProfile
from MBD_system.body.geometry.geometry_2D_profile import ContactGeometry2DProfile


class ContactRoughnessProfile(GeneralContact):
    """
    classdocs
    """

    def __init__(self, _type, body_id_i, body_id_j, u_iP=np.zeros(2), u_jP=np.zeros(2), properties_dict={}, parent=None):
        """
        
        :return:
        """
        super(ContactRoughnessProfile, self).__init__(_type, body_id_i, body_id_j, u_iP=u_iP, u_jP=u_jP, properties_dict=properties_dict, parent=parent)

        #   type of contact
        self._contact_type = "roughness profile"

        #   roughness profile list
        self.roughness_profile_list = []
        self.contact_geometry_list = []
        for body, body_id in zip(self.body_list, self.body_id_list):
            #   create roughness profile object
            roughness_profile = RoughnessProfile(body_id, body=body, _dict=properties_dict, parent=self)

            #   add roughness profile to list
            self.roughness_profile_list.append(roughness_profile)

            #   add geometry of surface roughness to list
            self.contact_geometry_list.append(roughness_profile.geometry.geom_data)

"""
Created on 18. mar. 2014

@author: lskrinjar (email: skrinjar.luka@gmail.com)
"""
import logging
import os
import time
from collections import OrderedDict
from pprint import pprint

import numpy as np


def create_list(filename, parent=None):
    """
    lists file joints.txt in folder - folder_path and returns the list of joints as objects
    Args:
        filename - absolute path to file
    Returns:
        contacts - list of contact objects
    """
    from MBD_system.contact.contact import Contact
    from MBD_system.contact.contact_general.contact_roughness_profile.contact_roughness_profile import ContactRoughnessProfile
    from MBD_system.contact.contact_pin_slot_clearance_joint.linear.pin_slot_clearance_joint_linear import PinSlotClearanceJointLinear
    from MBD_system.contact.contact_plane_sphere.contact_plane_sphere import ContactPlaneSphere
    from MBD_system.contact.contact_point_line.contact_point_line import ContactPointLine
    from MBD_system.contact.contact_revolute_clearance_joint.revolute_clearance_joint import RevoluteClearanceJoint
    from MBD_system.contact.contact_sphere_sphere.contact_sphere_sphere import ContactSphereSphere
    #    predefine empty dictionary
    dict_ = OrderedDict([])
    #    predefine contacts empty list
    contacts = []
    if os.path.isfile(filename):
        with open(filename, 'r') as file_:
        #    'r' - read
        #    'w' - write
        #    'a' - append
        #    'r+' - read and write
            for line in file_:
                line = line.replace('\n', '')
                __line = line.replace(' ', '')

                #   remove comment - after # sign
                if "#" in line:
                    line = line[:line.index("#")].strip()

                #   if line is empty or begins with #-comment, skip it
                if len(line.strip()) == 0 or line.startswith("#"):
                    pass
                else:
                    #    continue with reading of the file
                    if line.startswith("CONTACTS FILE-START") or line.startswith("#"):
                        pass

                    elif __line.startswith("type=") and len(dict_)>=4:
                        dict_["all_attributes_read"] = True

                    elif line.startswith("CONTACTS FILE-END"):
                        dict_["all_attributes_read"] = True

                    else:
                        key = line[0:line.index("=")].strip()

                        _val = line[line.index('=') + 1:].strip()

                        if _val.lower() == "true":
                            value = True
                        elif _val.lower() == "false":
                            value = False

                        else:
                            if key == "body_id_i" or key == "body_id_j":
                                value = int(line[line.index('=') + 1:].strip())

                            elif key == "u_iP" or key == "u_jP" or "u_" in key or "contact_area" in key or key == "n_i":
                                #   string to array
                                _array = line[line.index('=') + 1:].strip()
                                value = np.fromstring(_array, dtype=float, sep=',')

                            elif key == "K":
                                value = float(line[line.index('=') + 1:].strip())

                            else:
                                value = line[line.index('=') + 1:].strip()
                                try:
                                    #    string to integer
                                    value = int(line[line.index('=') + 1:].strip())
                                except:
                                    try:
                                        #    string to float
                                        value = float(line[line.index('=') + 1:].strip())
                                    except:
                                        value = line[line.index('=') + 1:].strip()
                        #    add dict item to dict with value
                        #   checks if key is already in dist it does not override current value, but adds/appends
                        #   new value to key item
                        if key in dict_:
                            #   check if key value is already a list
                            #   and if it is a list, then new value is appended to list
                            #   else old key value item is tranformed to list and
                            if isinstance(dict_[key], list):
                                dict_[key] = dict_[key].append(value)
                            else:
                                dict_[key] = [dict_[key], value]
                        else:
                            dict_[key] = value

                    if dict_.has_key("type") and dict_.has_key("all_attributes_read"):
                        if dict_["all_attributes_read"]:

                            #   general contact
                            if (dict_["type"]=="general") and (dict_.has_key("body_id_i")) and (dict_.has_key("body_id_j")):
                                #    construct contact object and append it to list
                                contact = Contact(dict_["type"], dict_["body_id_i"], dict_["body_id_j"], properties_dict=dict_, parent=parent)
                                contacts.append(contact)
                                #    the ordered dictionary in reinitialized to empty for properties of next contact to be read
                                dict_ = OrderedDict([])
                            #   revolute clearance joint
                            elif (dict_["type"] == "revolute clearance joint") and (dict_.has_key("body_id_i")) and (dict_.has_key("body_id_j")) and (dict_.has_key("u_iP")) and (dict_.has_key("u_jP")) and (dict_.has_key("R0_i")) and (dict_.has_key("R0_j")):
                                #    construct revolute clearance joint (contact type) and append it to list
                                contact = RevoluteClearanceJoint(dict_["type"], dict_["body_id_i"], dict_["body_id_j"], dict_["u_iP"], dict_["u_jP"], dict_["R0_i"], dict_["R0_j"], properties_dict=dict_, parent=parent)
                                contacts.append(contact)
                                #    the ordered dictionary in reinitialized to empty for properties of next contact to be read
                                dict_ = OrderedDict([])
                            #   sphere-sphere contact
                            elif (dict_["type"] == "contact sphere-sphere") and (dict_.has_key("body_id_i")) and (dict_.has_key("body_id_j")) and (dict_.has_key("R_i")) and (dict_.has_key("R_j")):
                                contact = ContactSphereSphere(dict_["type"], dict_["body_id_i"], dict_["body_id_j"], dict_["R_i"], dict_["R_j"], properties_dict=dict_, parent=parent)
                                contacts.append(contact)
                                #    the ordered dictionary in reinitialized to empty for properties of next contact to be read
                                dict_ = OrderedDict([])
                            #   plane-sphere contact
                            elif (dict_["type"] == "contact plane-sphere") and (dict_.has_key("body_id_i")) and (dict_.has_key("body_id_j")) and (dict_.has_key("R0_j")) and (dict_.has_key("n_i")):
                                contact = ContactPlaneSphere(dict_["type"], dict_["body_id_i"], dict_["body_id_j"], u_iP=dict_["u_iP"], R0_j=dict_["R0_j"], n_i=dict_["n_i"], properties_dict=dict_, parent=parent)
                                contacts.append(contact)
                                #    the ordered dictionary in reinitialized to empty for properties of next contact to be read
                                dict_ = OrderedDict([])
                            #   pin-slot clearance joint linear
                            elif (dict_["type"]=="pin-slot clearane joint linear") and (dict_.has_key("body_id_i")) and (dict_.has_key("body_id_j")) and (dict_.has_key("u_iP")) and (dict_.has_key("u_iR")) and (dict_.has_key("u_jP")) and (dict_.has_key("h0_iP")) and (dict_.has_key("R0_j")):
                                contact = PinSlotClearanceJointLinear(dict_["type"], dict_["body_id_i"], dict_["body_id_j"], dict_["u_iP"], dict_["u_iR"], dict_["u_jP"], dict_["h0_iP"], dict_["R0_j"], properties_dict=dict_, parent=parent)
                                #    the ordered dictionary in reinitialized to empty for properties of next contact to be read
                                contacts.append(contact)
                                dict_ = OrderedDict([])
                            #   pin-slot clearance joint radial
                            elif (dict_["type"]=="contact plane-radial"):
                                print dict_["type"] + "Under construction!"
                            #   prismatic clearance joint radial
                            elif (dict_["type"]=="prismatic clearance joint"):
                                print dict_["type"] + "Under construction!"
                            #   point-line contact
                            elif (dict_["type"] == "contact point-line") and (dict_.has_key("body_id_i")) and (dict_.has_key("body_id_j")):
                                contact = ContactPointLine(dict_["type"], dict_["body_id_i"], dict_["body_id_j"], dict_["u_iP"], dict_["u_jP"], dict_["u_jR"], properties_dict=dict_, parent=parent)
                                #    the ordered dictionary in reinitialized to empty for properties of next contact to be read
                                contacts.append(contact)
                                dict_ = OrderedDict([])
                            elif (dict_["type"] == "roughness profile") and (dict_.has_key("body_id_i")) and (dict_.has_key("body_id_j")):
                                contact = ContactRoughnessProfile(dict_["type"], dict_["body_id_i"], dict_["body_id_j"], dict_["u_iP"], dict_["u_jP"], properties_dict=dict_, parent=parent)
                                contacts.append(contact)
                                #    the ordered dictionary in reinitialized to empty for properties of next contact to be read
                                dict_ = OrderedDict([])
                            else:
                                print dict_
                                print "type =", dict_["type"]
                                print "Contact not created!"
                                # raise ValueError,

                    if __line.startswith("type="):
                        key = line[0:line.index("=")].strip()
                        value = line[line.index("=")+1:].strip()
                        dict_[key] = value

                    if line.startswith("CONTACTS FILE-END"):
                        file_.close()
                        break
        logging.getLogger("DyS_logger").info("File %s found! Contacts created successfully!", filename)
    else:
        logging.getLogger("DyS_logger").info("File %s not found! No contacts created!", filename)

    return contacts

        
if __name__ == "__main__":
    t1 = time.time()
    contact_list = create_list("contacts.dat")
    print "-------------------------"
    print "number of contacts =", len(contact_list)
#     print time.time() - t1
    for contact in contact_list:
        pprint(vars(contact))

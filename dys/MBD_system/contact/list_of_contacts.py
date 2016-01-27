'''
Created on 18. mar. 2014

@author: lskrinjar (email: skrinjar.luka@gmail.com)
'''
import logging
import os
import time
from collections import OrderedDict
from pprint import pprint

import numpy as np

from MBD_system.contact.contact import Contact
from MBD_system.contact.contact_revolute_clearance_joint.revolute_clearance_joint import RevoluteClearanceJoint


def create_list(filename, parent=None):
    """
    lists file joints.txt in folder - folder_path and returns the list of joints as objects
    Args:
        filename - absolute path to file
    Returns:
        contacts - list of contact objects
    """
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

                        if key == "body_i" or key == "body_j":
                            value = int(line[line.index('=') + 1:].strip())

                        elif key == "u_iP" or key == "u_jP" or "contact_area" in key:
                            #   string to array
                            _array = line[line.index('=') + 1:].strip()
                            value = np.fromstring(_array, dtype=float, sep=',')
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
                        dict_[key] = value

                    if dict_.has_key("type") and dict_.has_key("all_attributes_read"):
                        if dict_["all_attributes_read"]:
                            if (dict_.has_key("body_i")) and (dict_.has_key("body_j")) and (dict_["type"]=="general"):
                                #    construct contact object and append it to list
                                contact = Contact(dict_["body_i"], dict_["body_j"], dict_["type"], properties_dict=dict_, parent=parent)
                                contacts.append(contact)
                                #    the ordered dictionary in reinitialized to empty for properties of next contact to be read
                                dict_ = OrderedDict([])
                            if (dict_.has_key("body_i")) and (dict_.has_key("body_j")) and (dict_.has_key("u_iP")) and (dict_.has_key("u_jP")) and (dict_.has_key("R_i")) and (dict_.has_key("R_j")) and (dict_["type"]=="revolute clearance joint"):
                                #    construct revolute clearance joint (contact type) and append it to list
                                contact = RevoluteClearanceJoint(dict_["body_i"], dict_["body_j"], dict_["u_iP"], dict_["u_jP"], dict_["R_i"], dict_["R_j"], dict_["type"], properties_dict=dict_, parent=parent)
                                contacts.append(contact)
                                #    the ordered dictionary in reinitialized to empty for properties of next contact to be read
                                dict_ = OrderedDict([])

                    if __line.startswith("type="):
                        key = line[0:line.index("=")].strip()
                        value = line[line.index("=")+1:].strip()
                        dict_[key] = value

                    if line.startswith("CONTACTS FILE-END"):
                        file_.closed
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

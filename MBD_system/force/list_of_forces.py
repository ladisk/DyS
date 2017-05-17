"""
Created on 18. mar. 2014

@author: lskrinjar (email: skrinjar.luka@gmail.com)
"""
import logging
import os
import re
from collections import OrderedDict
from pprint import pprint

import numpy as np


def create_list(filename, parent=None):
    """
    Lists file joints.txt in folder - folder_path and returns the list of joints as objects
    :param filename:    absolute path to file
    :returns forces:    list of force objects
    """
    from force import Force
    #    predefine empty dictionary
    dict_ = OrderedDict([])
    #   predefine contacts empty list
    forces = []
    if os.path.isfile(filename):
        lines = open(filename).read().split('\n')
#         with open(filename, "r") as lines:
        #    'r' - read
        #    'w' - write
        #    'a' - append
        #    'r+' - read and write
#         content = lines.readlines()
        #    close file
        for line in lines:

            if "#" in line:
                line = line[:line.index("#")].strip()

            if line.startswith("FORCES FILE-START") or line.startswith("#"):
                continue

            if line.startswith("force_name") or line.startswith("FORCES FILE-END"):
                if len(dict_) >= 4:
                    if "node_id" not in dict_:
                        dict_["node_id"] = None

                    if "u_iP_f" not in dict_:
                        dict_["u_iP_f"] = np.zeros(2)

                    force = Force(dict_["body_id"], force_name=dict_["force_name"], Fx=dict_["Fx"], Fy=dict_["Fy"], Mz=dict_["Mz"], u_iP_f=dict_["u_iP_f"], node_id=dict_["node_id"], dict=dict_, parent=parent)
                    forces.append(force)
                    #   cleanup dictionary
                    dict_ = OrderedDict([])

            if line == "FORCES FILE-END":
                break

            else:
                match_equal = re.search(r"=", line)
                key = line[0:match_equal.start()].strip()
                match_newline = re.search(r"\n", line)
                val = line[match_equal.start()+1:match_newline].strip()

                if key == "body_id":
                    val = int(val)

                elif key == "u_iP_f":
                    val = np.array(val.split(','), dtype="float64")

                elif key == "node_id":
                    val = int(val)

                elif key == "Fx" or key == "Fy" or key == "Mz":
                    try:
                        val = float(val)
                    except:
                        pass
                else:
                    try:
                        val = float(val)
                    except:
                        pass

                if val == "True":
                    val = True

                if val == "False":
                    val = False

                #   add attributes to dictionary
                dict_[key] = val

        # logging.getLogger("DyS_logger").info("File %s found! Forces created successfully!", filename)
        
    else:
        msg = "File %s not found! No forces created!", filename
        print msg
        logging.getLogger("DyS_logger").info(msg)
        
    return forces

        
if __name__ == "__main__":
    list = create_list("forces.txt")
    for object in list:
        pprint(vars(object))


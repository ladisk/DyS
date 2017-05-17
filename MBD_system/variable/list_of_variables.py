"""
Created on 18. mar. 2014

@author: lskrinjar (email: skrinjar.luka@gmail.com)
"""
import logging
import os
import re
from collections import OrderedDict
import time
from pprint import pprint

import numpy as np


def create_list(filename, parent=None):
    """
    lists file joints.txt in folder - folder_path and returns the list of joints as objects
    in:
        filename - absolute path to file
    returns:
        springs - list of contact objects
    """
    from variable import Variable

    #    predefine empty dictionary
    dict_ = OrderedDict([])

    #   predefine empty list
    variables = []

    if os.path.isfile(filename):
        with open(filename, 'r') as file_:
        #    'r' - read
        #    'w' - write
        #    'a' - append
        #    'r+' - read and write
            
            #    joint type - predefine params
            spring_name_ = None
            spring_type_ = None
            k_ = 0
            c_ = 0
            k_t_ = 0
            c_t_ = 0
            body_i_ = None
            body_j_ = None
            u_iP_ = None
            u_jP_ = None
            l_0_ = None
            F0_ = None
            
            for line in file_:
                if line.startswith("VARIABLES-START") or line.startswith("#"):
                    continue

                if line.startswith("mbd_item_type") or line.startswith("VARIABLES-END"):
                    if len(dict_) == 6:
                        variable = Variable(MBD_item_type=dict_["mbd_item_type"],
                                            item_id=dict_["item_id"],
                                            variable=dict_["variable"],
                                            x_n=dict_["x_n"],
                                            x_l=dict_["x_l"],
                                            x_u=dict_["x_u"],
                                            parent=parent)
                        variables.append(variable)

                        dict_ = OrderedDict([])

                if line.startswith("VARIABLES-END"):
                    break

                else:
                    match_equal = re.search(r"=", line)
                    key = line[0:match_equal.start()].strip()
                    match_newline = re.search(r"\n", line)
                    val = line[match_equal.start()+1:match_newline.start()].strip()

                if key == "mbd_item_type":
                    val = str(val)

                elif key == "item_id":
                    val = int(val)

                elif key == "variable":
                    val = str(val)

                elif key in ["x_n", "x_u", "x_l"]:
                    if "," in val:
                        val = np.array(val.split(','), dtype="float32")
                    else:
                        val = float(val)

                #   add attributes to dictionary
                dict_[key] = val

            # logging.getLogger("DyS_logger").info("File %s found! Springs created successfully!", filename)
            file_.close()
    else:
        pass
        # logging.getLogger("DyS_logger").info("File %s not found! No variables created!", filename)
                
    return variables

        
if __name__ == "__main__":
    t1 = time.time()
    spring_list = create_list("springs.txt")
    print time.time() - t1
    for spring in spring_list:
        pprint(vars(spring))
    

    
    

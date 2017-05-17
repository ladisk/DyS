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
from sympy.polys.tests.test_densearith import F_0


def create_list(filename, parent=None):
    """
    lists file joints.txt in folder - folder_path and returns the list of joints as objects
    in:
        filename - absolute path to file
    returns:
        springs - list of contact objects
    """
    from spring_translational import SpringTranslational
    from spring_rotational import SpringRotational
    from spring_translational_on_flexible_body_to_ground import SpringTranslationalOnFlexibleBodyToGround
    #    predefine empty dictionary
    dict_ = OrderedDict([])

    #   predefine empty list
    springs = []

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
            ksi_ = None
            element_id_ = None
            element_ksi_ = None
            joint_id_ = None
            
            for line in file_:
                if line.startswith("SPRINGS FILE-START") or line.startswith("#"):
                    continue

                if line.startswith("spring_name") or line.startswith("SPRINGS FILE-END"):
                    spring = None

                    if len(dict_) >= 7 and dict_["spring_type"] == "translational":
                        spring = SpringTranslational(dict_["spring_name"], dict_["body_id_i"], dict_["body_id_j"], dict_["u_iP"], dict_["u_jP"], k=dict_["k"], c=dict_["c"], l_0=dict_["l_0"], properties_dict=dict_, parent=parent)

                    elif (len(dict_) >= 5 and dict_["spring_type"] == "rotational") or ("joint_id" in dict_):
                        spring = SpringRotational(dict_["spring_name"], dict_["body_id_i"], dict_["body_id_j"], k_t=dict_["k_t"], c_t=dict_["c_t"], properties_dict=dict_, parent=parent)

                    elif len(dict_) >= 9 and dict_["spring_type"] == "translational on flexible body to ground":
                        spring = SpringTranslationalOnFlexibleBodyToGround(dict_["spring_name"], dict_["body_id_i"], u_jP_CAD=dict_["u_jP"], k=dict_["k"], c=dict_["c"], l_0=dict_["l_0"], properties_dict=dict_, parent=parent)

                if spring is not None:
                    springs.append(spring)
                    spring = None

                if line.startswith("SPRINGS FILE-END"):
                    break

                else:
                    match_equal = re.search(r"=", line)
                    key = line[0:match_equal.start()].strip()
                    match_newline = re.search(r"\n", line)
                    val = line[match_equal.start()+1:match_newline.start()].strip()

                if key in ["body_id_i", "body_id_j", "joint_id"]:
                    try:
                        val = int(val)
                    except:
                        pass

                elif key in ["u_iP", "u_jP", "I_r_ij_0", "color"]:
                    val = np.array(val.split(','), dtype="float64")

                elif key in ["k", "c", "k_t", "c_t", "l_0"]:
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


            # logging.getLogger("DyS_logger").info("File %s found! Springs created successfully!", filename)
        file_.close()
    else:
        logging.getLogger("DyS_logger").info("File %s not found! No springs created!", filename)
                
    return springs

        
if __name__ == "__main__":
    t1 = time.time()
    spring_list = create_list("springs.txt")
    print time.time() - t1
    for spring in spring_list:
        pprint(vars(spring))
    

    
    

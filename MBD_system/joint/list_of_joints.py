"""
Created on 18. mar. 2014

@author: lskrinjar (email: skrinjar.luka@gmail.com)
"""
import logging
import os
from pprint import pprint
import re
import time


from joint import Joint
import numpy as np


def create_list(joints_=None, filename=None, parent=None):
    """
    lists file joints.txt in folder - folder_path and returns the list of joints as objects
    :param
        folder_path - absolute path of folder
        folder_name - name of the folder
    """

#     os.chdir(os.path.abspath(__file__))
    _parent = parent

    if os.path.isfile(filename):
        with open(filename, 'r') as file_:
        #    'r' - read
        #    'w' - write
        #    'a' - append
        #    'r+' - read and write
            
            #    joint type - predefine params
            joint_type_ = None
            body_i_ = None
            body_j_ = None
            u_iP_ = None
            u_jP_ = None
            u_iQ_ = None
            
            for line in file_:
                if line.startswith("JOINTS FILE-START"):
                    continue
                elif line.startswith("JOINTS FILE-END"):
                    file_.close()
                    break
                
                elif line.startswith("joint_type="):
                    match_equal = re.search(r"=", line)
                    match_newline = re.search(r"\n", line)
                    joint_type_ = line[match_equal.start() + 1:match_newline.start()]
    
                elif line.startswith("body_i"):
                    match_equal = re.search(r"=", line)
                    match_newline = re.search(r"\n", line)
                    val_ = line[match_equal.start() + 2:match_newline.start()]
                    try:
                        body_i_ = np.int(val_)
    #                         print "int =", body_i_
                    except:
                        body_i_ = str(val_)
    #                         print "str =", body_i_    
    
                elif line.startswith("body_j"):
                    match_equal = re.search(r"=", line)
                    match_newline = re.search(r"\n", line)
                    val_ = line[match_equal.start() + 2:match_newline.start()]
                    try:
                        body_j_ = np.int(val_)
    #                         print "int =", body_i_
                    except:
                        body_j_ = str(val_)
    #                         print "str =", body_i_  
                    
                elif line.startswith("u_iP"):
                    match_equal = re.search(r"=", line)
                    match_newline = re.search(r"\n", line)
                    u_iP_ = np.array(line[match_equal.start() + 2:match_newline.start()].split(','), dtype="float64")
                    
                elif line.startswith("u_iQ"):
                    match_equal = re.search(r"=", line)
                    match_newline = re.search(r"\n", line)
                    u_iQ_ = np.array(line[match_equal.start() + 2:match_newline.start()].split(','), dtype="float64")
                    
                    
                elif line.startswith("u_jP"):
                    match_equal = re.search(r"=", line)
                    match_newline = re.search(r"\n", line)
                    u_jP_ = np.array(line[match_equal.start() + 2:match_newline.start()].split(','), dtype="float64")
    
                #   create revolute joint
                if (joint_type_ != None) and (body_i_ != None) and (body_j_ != None) and (u_iP_ != None) and (u_jP_ != None):
                    #   create joint object
                    joint_ = Joint(joint_type_, body_i_, body_j_, u_iP_CAD=u_iP_, u_jP_CAD=u_jP_, parent=_parent)
                    joints_.append(joint_)
    
                    #    joint type
                    joint_type_ = None
                    body_i_ = None
                    body_j_ = None
                    u_iP_ = None
                    u_jP_ = None

                #   create prismatic joint
                elif (joint_type_ != None) and (body_i_ != None) and (body_j_ != None) and (u_iP_ != None) and (u_jP_ != None) and (u_iQ_ != None):
                    #    create joint object here
                    joint_ = Joint(joint_type_, body_i_, body_j_, u_iP_CAD=u_iP_, u_jP_CAD=u_jP_, u_iQ=u_iQ_, parent=_parent)
                    joints_.append(joint_)
                    
                    #    joint type
                    joint_type_ = None
                    body_i_ = None
                    body_j_ = None
                    u_iP_ = None
                    u_jP_ = None
        # logging.getLogger("DyS_logger").info("File %s found! Joints created successfully!", filename)
    
    else:
        logging.getLogger("DyS_logger").info("File %s not found! No joints created!", filename)

    return joints_

        
if __name__ == "__main__":
    t1 = time.time()
    joint_list = create_list("joints.txt")
    print time.time() - t1
    print "joint_list =", joint_list
    for joint in joint_list:
        pprint(vars(joint))
    

    
    

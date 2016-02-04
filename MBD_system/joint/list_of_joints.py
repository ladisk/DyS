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
from joint_revolute import JointRevolute
from joint_prismatic import JointPrismatic
from joint_fixed import JointFixed
import numpy as np


def create_list(joints, filename, parent=None):
    """
    Create list of joints from file joint.dat
    :param filename:    absolute path to file
    :param parent:      parent to joint object
    """
    if os.path.isfile(filename):
        with open(filename, 'r') as file_:
        #    'r' - read
        #    'w' - write
        #    'a' - append
        #    'r+' - read and write
            
            #    joint type - predefine params
            _joint_type = None
            _body_i = None
            _body_j = None
            _u_iP = None
            _u_jP = None
            _u_iQ = None
            
            for line in file_:
                if line.startswith("JOINTS FILE-START"):
                    continue
                elif line.startswith("JOINTS FILE-END"):
                    file_.close()
                    break
                
                elif line.startswith("joint_type="):
                    match_equal = re.search(r"=", line)
                    match_newline = re.search(r"\n", line)
                    _joint_type = line[match_equal.start() + 1:match_newline.start()]
    
                elif line.startswith("body_i"):
                    match_equal = re.search(r"=", line)
                    match_newline = re.search(r"\n", line)
                    val_ = line[match_equal.start() + 2:match_newline.start()]
                    try:
                        _body_i = np.int(val_)
    #                         print "int =", _body_i
                    except:
                        _body_i = str(val_)
    #                         print "str =", _body_i    
    
                elif line.startswith("body_j"):
                    match_equal = re.search(r"=", line)
                    match_newline = re.search(r"\n", line)
                    val_ = line[match_equal.start() + 2:match_newline.start()]
                    try:
                        _body_j = np.int(val_)
    #                         print "int =", _body_i
                    except:
                        _body_j = str(val_)
    #                         print "str =", _body_i  
                    
                elif line.startswith("u_iP"):
                    match_equal = re.search(r"=", line)
                    match_newline = re.search(r"\n", line)
                    _u_iP = np.array(line[match_equal.start() + 2:match_newline.start()].split(','), dtype="float64")
                    
                elif line.startswith("u_iQ"):
                    match_equal = re.search(r"=", line)
                    match_newline = re.search(r"\n", line)
                    _u_iQ = np.array(line[match_equal.start() + 2:match_newline.start()].split(','), dtype="float64")

                elif line.startswith("u_jP"):
                    match_equal = re.search(r"=", line)
                    match_newline = re.search(r"\n", line)
                    _u_jP = np.array(line[match_equal.start() + 2:match_newline.start()].split(','), dtype="float64")
    
                #   create revolute joint
                if (_joint_type == "revolute") and (_body_i != None) and (_body_j != None) and (_u_iP != None) and (_u_jP != None):
                    #   create joint object
                    joint = JointRevolute(_joint_type, _body_i, _body_j, u_iP_CAD=_u_iP, u_jP_CAD=_u_jP, parent=parent)
                    joints.append(joint)
    
                    #    reset variables
                    _joint_type = _body_i = _body_j = _u_iP = _u_jP = None

                #   create prismatic joint
                elif (_joint_type == "prismatic") and (_body_i != None) and (_body_j != None) and (_u_iP != None) and (_u_jP != None) and (_u_iQ != None):
                    #    create joint object
                    print "input"
                    print "_u_iP =", _u_iP
                    print "_u_jP =", _u_jP
                    print "_u_iQ =", _u_iQ
                    joint = JointPrismatic(_joint_type, _body_i, _body_j, u_iP_CAD=_u_iP, u_jP_CAD=_u_jP, u_iQ_CAD=_u_iQ, parent=parent)
                    joints.append(joint)
                    
                    #    reset variables
                    _joint_type = _body_i = _body_j = _u_iP = _u_jP = _u_iQ = None

                #   create fixed joint
                elif (_joint_type == "fixed") and (_body_i != None) and (_body_j != None) and (_u_iP != None) and (_u_jP != None):
                    #   create joint object here
                    joint = JointFixed(_joint_type, _body_i, _body_j, u_iP_CAD=_u_iP, u_jP_CAD=_u_jP, parent=parent)
                    joints.append(joint)

                    #    reset variables
                    _joint_type = _body_i = _body_j = _u_iP = _u_jP = None

        # logging.getLogger("DyS_logger").info("File %s found! Joints created successfully!", filename)
    
    else:
        logging.getLogger("DyS_logger").info("File %s not found! No joints created!", filename)

    return joints

        
if __name__ == "__main__":
    t1 = time.time()
    joint_list = create_list("joints.txt")
    print time.time() - t1
    print "joint_list =", joint_list
    for joint in joint_list:
        pprint(vars(joint))
    

    
    

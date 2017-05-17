"""
Created on 18. mar. 2014

@author: lskrinjar (email: skrinjar.luka@gmail.com)
"""
import logging
import os
import re
import time
from collections import OrderedDict
from pprint import pprint

import numpy as np


def create_list(joints, filename, parent=None):
    """
    Create list of joints from file joint.dat
    :param filename:    absolute path to file
    :param parent:      parent to joint object
    """
    from joint_fixed import JointFixed
    from joint_prismatic import JointPrismatic
    from joint_revolute import JointRevolute
    from support_fixed import SupportFixed
    from support_hinged import SupportHinged
    from support_roller import SupportRoller
    from joint_revolute_rigid_flexible import JointRevoluteRigidFlexible
    from joint_rigid_rigid_flexible import JointRigidRigidFlexible
    from joint_revolute_flexible_flexible import JointRevoluteFlexibleFlexible
    from joint_rigid_flexible_flexible import JointRigidFlexibleFlexible
    from joint_rigid_point_mass_flexible import JointRigidPointMassFlexible
    from joint_rigid_point_mass_flexible_with_rotation import JointRigidPointMassFlexibleWithRotation
    from joint_fixed_point_mass_rigid import JointFixedPointMassRigid

    if os.path.isfile(filename):
        with open(filename, 'r') as file_:
            #   'r' - read
            #   'w' - write
            #   'a' - append
            #   'r+' - read and write
            
            #    joint type - predefine params
            _joint_type = None
            _body_i = None
            _body_j = None
            _u_iP = None
            _u_jP = None
            _u_iQ = None
            _node_id_i = None
            _node_id_j = None
            _u_iR = None
            #    predefine empty dictionary
            dict_ = OrderedDict([])
            
            for line in file_:
                line = line.rstrip()

                #   remove comment - after # sign
                if "#" in line:
                    line = line[:line.index("#")]

                # print line
                if line.startswith("JOINTS FILE-START") or line.startswith("#") or line.startswith("JOINTS FILE-END"):
                    pass

                elif line.startswith("body_id_i"):
                    match_equal = re.search(r"=", line)
#                     match_newline = re.search(r"\n", line)
                    val_ = line[match_equal.start() + 2:]
                    
                    try:
                        _body_i = np.int(val_)
    #                         print "int =", _body_i
                    except:
                        _body_i = str(val_)
    #                         print "str =", _body_i

                elif line.startswith("body_id_j"):
                    match_equal = re.search(r"=", line)
#                     match_newline = re.search(r"\n", line)
                    val_ = line[match_equal.start() + 2:]
                    try:
                        _body_j = np.int(val_)
    #                         print "int =", _body_i
                    except:
                        _body_j = str(val_)
    #                         print "str =", _body_i

                elif line.startswith("u_iP"):
                    match_equal = re.search(r"=", line)
                    _u_iP = np.array(line[match_equal.start() + 2:].split(','), dtype="float64")#match_newline.start()].split(',')

                elif line.startswith("u_iQ"):
                    match_equal = re.search(r"=", line)
                    _u_iQ = np.array(line[match_equal.start() + 2:].split(','), dtype="float64")#match_newline.start()].split(',')

                elif line.startswith("u_jP"):
                    match_equal = re.search(r"=", line)
                    _u_jP = np.array(line[match_equal.start() + 2:].split(','), dtype="float64")

                elif line.startswith("u_iR"):
                    match_equal = re.search(r"=", line)
                    _u_iR = np.array(line[match_equal.start() + 2:].split(','), dtype="float64")

                elif line.startswith("node_id_jk"):
                    match_equal = re.search(r"=", line)
                    _node_id_j = int(line[match_equal.start() + 2:])

                elif line.startswith("node_id_i"):
                    match_equal = re.search(r"=", line)
                    _node_id_i = int(line[match_equal.start() + 2:])

                elif line.startswith("node_id_j"):
                    match_equal = re.search(r"=", line)
                    _node_id_j = int(line[match_equal.start() + 2:])

                else:
                    #    key
                    key = line[0:line.index("=")]
                    if key == "name":
                        key = "_" + key

                    #    value
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

                    if key != "joint_type":
                        #   add dict item to dict with value
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

                if line.startswith("joint_type") or line.startswith("JOINTS FILE-END") or line.startswith("#"):
                    #   create revolute joint
                    if (_joint_type == "revolute") and (_body_i is not None) and (_body_j is not None) and (_u_iP is not None) and (_u_jP is not None):
                        #   create joint object
                        joint = JointRevolute(_joint_type, _body_i, _body_j, u_iP_CAD=_u_iP, u_jP_CAD=_u_jP, properties_dict=dict_, parent=parent)
                        joints.append(joint)

                        #    reset variables
                        _joint_type = _body_i = _body_j = _u_iP = _u_jP = None
                        dict_ = OrderedDict([])

                    #   create prismatic joint
                    elif (_joint_type == "prismatic") and (_body_i is not None) and (_body_j is not None) and (_u_iP is not None) and (_u_jP is not None) and (_u_iQ is not None):
                        #    create joint object
                        joint = JointPrismatic(_joint_type, _body_i, _body_j, u_iP_CAD=_u_iP, u_jP_CAD=_u_jP, u_iQ_CAD=_u_iQ, properties_dict=dict_, parent=parent)
                        joints.append(joint)

                        #    reset variables
                        _joint_type = _body_i = _body_j = _u_iP = _u_jP = _u_iQ = None
                        dict_ = OrderedDict([])

                    #   create fixed joint
                    elif (_joint_type == "fixed") and (_body_i is not None) and (_body_j is not None) and (_u_iP is not None) and (_u_jP is not None):
                        #   create joint object here
                        joint = JointFixed(_joint_type, _body_i, _body_j, u_iP_CAD=_u_iP, u_jP_CAD=_u_jP, properties_dict=dict_, parent=parent)
                        joints.append(joint)

                        #    reset variables
                        _joint_type = _body_i = _body_j = _u_iP = _u_jP = None
                        dict_ = OrderedDict([])

                    elif (_joint_type == "fixed support") and (_body_i is not None) and (_body_j is not None) and (_node_id_i is not None) and (_node_id_j is not None):
                        #   create fixed support object for ANCF body
                        joint = SupportFixed(_body_i, _body_j, node_id_i=_node_id_i, node_id_j=_node_id_j, properties_dict=dict_, parent=parent)
                        joints.append(joint)

                        #    reset variables
                        _joint_type = _body_i = _body_j = _node_id_i = _node_id_j = None
                        dict_ = OrderedDict([])

                    elif (_joint_type == "hinged support") and (_body_i is not None) and (_body_j is not None) and (_node_id_i is not None) and (_node_id_j is not None):
                        #   create fixed support object for ANCF body
                        joint = SupportHinged(_body_i, _body_j, node_id_i=_node_id_i, node_id_j=_node_id_j, properties_dict=dict_, parent=parent)
                        joints.append(joint)

                        #    reset variables
                        _joint_type = _body_i = _body_j = _node_id_i = _node_id_j = None
                        dict_ = OrderedDict([])

                    elif (_joint_type == "revolute joint rigid-flexible") and (_body_i is not None) and (_body_j is not None) and (_u_iP is not None) and (_node_id_j is not None):
                        #   create revolute joint object between rigid and flexible body
                        joint = JointRevoluteRigidFlexible(_body_i, _body_j, u_iP=_u_iP, node_id_j=_node_id_j, properties_dict=dict_, parent=parent)
                        joints.append(joint)

                        #    reset variables
                        _joint_type = _body_i = _body_j = _node_id_i = _u_iP = _node_id_j = None
                        dict_ = OrderedDict([])

                    elif (_joint_type == "rigid joint rigid-flexible") and (_body_i is not None) and (_body_j is not None) and (_u_iP is not None) and (_node_id_j is not None) and (_u_iR is not None):
                        #   create revolute joint object between rigid and flexible body
                        joint = JointRigidRigidFlexible(_body_i, _body_j, node_id_j=_node_id_j, u_iP=_u_iP, u_iR=_u_iR, properties_dict=dict_, parent=parent)
                        joints.append(joint)

                        #    reset variables
                        _joint_type = _body_i = _body_j = _u_iP = _node_id_j = None
                        dict_ = OrderedDict([])

                    elif (_joint_type == "revolute joint flexible-flexible") and (_body_i is not None) and (_body_j is not None) and (_node_id_i is not None) and (_node_id_j is not None):
                        # create revolute joint object between rigid and flexible body
                        joint = JointRevoluteFlexibleFlexible(_body_i, _body_j, node_id_i=_node_id_i, node_id_j=_node_id_j, properties_dict=dict_, parent=parent)
                        joints.append(joint)

                        #    reset variables
                        _joint_type = _body_i = _body_j = _node_id_i = _node_id_j = None
                        dict_ = OrderedDict([])

                    elif (_joint_type == "rigid joint flexible-flexible") and (_body_i is not None) and (_body_j is not None) and (_node_id_i is not None) and (_node_id_j is not None):
                        joint = JointRigidFlexibleFlexible(_body_i, _body_j, node_id_i=_node_id_i, node_id_j=_node_id_j, properties_dict=dict_, parent=parent)
                        joints.append(joint)

                        #    reset variables
                        _joint_type = _body_i = _body_j = _node_id_i = _node_id_j = None
                        dict_ = OrderedDict([])

                    elif (_joint_type == "rigid joint point mass-flexible") and (_body_i is not None) and (_body_j is not None) and (_u_iP is not None) and (_node_id_j is not None):
                        joint = JointRigidPointMassFlexible(_body_i, _body_j, u_iP=_u_iP, node_id_j=_node_id_j, properties_dict=dict_, parent=parent)
                        joints.append(joint)

                        #    reset variables
                        _joint_type = _body_i = _body_j = _u_iP = _node_id_j = None
                        dict_ = OrderedDict([])

                    elif (_joint_type == "rigid joint point mass-flexible with rotation") and (_body_i is not None) and (_body_j is not None) and (_u_iP is not None) and (_node_id_j is not None):
                        joint = JointRigidPointMassFlexibleWithRotation(_body_i, _body_j, u_iP=_u_iP, node_id_j=_node_id_j, properties_dict=dict_, parent=parent)
                        joints.append(joint)

                        #    reset variables
                        _joint_type = _body_i = _body_j = _u_iP = _node_id_j = None
                        dict_ = OrderedDict([])

                    elif (_joint_type == "fixed joint point mass-rigid") and (_body_i is not None) and (_body_j is not None) and (_u_iP is not None) and (_u_jP is not None):
                        joint = JointFixedPointMassRigid(_body_i, _body_j, u_iP=_u_iP, u_jP=_u_jP, properties_dict=dict_, parent=parent)
                        joints.append(joint)

                        #    reset variables
                        _joint_type = _body_i = _body_j = _u_iP = _u_jP = None
                        dict_ = OrderedDict([])


                    elif _joint_type is None and _body_i is None and _body_j is None and _u_iP is None and _node_id_j is None:
                        pass

                    else:
                        print "Joint construction error! Not all joint parameters defined!"
                        print "_joint_type =", _joint_type
                        print "_body_i =", _body_i
                        print "_body_j =", _body_j

                if line.startswith("joint_type"):
                    line = line.replace(" = ", "=").replace("= ", "=").replace(" =", "=")
                    match_equal = re.search(r"=", line)
                    _joint_type = line.strip()[match_equal.start() + 1:]

                if line.startswith("JOINTS FILE-END"):
                    file_.close()
                    break

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
    

    
    

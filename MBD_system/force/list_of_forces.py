"""
Created on 18. mar. 2014

@author: lskrinjar (email: skrinjar.luka@gmail.com)
"""
import os
import re
import logging
from pprint import pprint
import time
import numpy as np


from force import Force


def create_list(filename, parent=None):
    """
    Lists file joints.txt in folder - folder_path and returns the list of joints as objects
    in:
        filename - absolute path to file
    returns:
        forces - list of contact objects
    """
    forces = []

    try:
        with open(filename, "r") as lines:
        #    'r' - read
        #    'w' - write
        #    'a' - append
        #    'r+' - read and write
            
            #    force - predefine params
            force_name_ = None
            body_id_ = None
            Fx_ = None
            Fy_ = None
            Mz_ = None
            u_iP_f_ = None
            
            for line in lines:
                if line.startswith("FORCES FILE-START"):
                    continue
                
                elif line.startswith("FORCES FILE-END"):
                    filename.close
                    break
                
                elif line.startswith("force_name="):
                    match_equal = re.search(r"=", line)
                    match_newline = re.search(r"\n", line)
                    force_name_ = line[match_equal.start() + 1:match_newline.start()]
                    
                elif line.startswith("body_id"):
                    match_equal = re.search(r"=", line)
                    match_newline = re.search(r"\n", line)
                    body_id_ = np.int(line[match_equal.start() + 2:match_newline.start()])
    
                elif line.startswith("Fx"):
                    match_equal = re.search(r"=", line)
                    match_newline = re.search(r"\n", line)
                    try:
                        Fx_ = np.float(line[match_equal.start() + 2:match_newline.start()])
                    except:
                        Fx_ = line[match_equal.start() + 2:match_newline.start()]
                    
                elif line.startswith("Fy"):
                    match_equal = re.search(r"=", line)
                    match_newline = re.search(r"\n", line)
                    try:
                        Fy_ = np.float(line[match_equal.start() + 2:match_newline.start()])
                    except:
                        Fy_ = line[match_equal.start() + 2:match_newline.start()]
                
                elif line.startswith("Mz"):
                    match_equal = re.search(r"=", line)
                    match_newline = re.search(r"\n", line)
                    try:
                        Mz_ = np.float(line[match_equal.start() + 2:match_newline.start()])
                    except:
                        Fy_ = line[match_equal.start() + 2:match_newline.start()]
                    
    
                elif line.startswith("u_iP_f"):
                    match_equal = re.search(r"=", line)
                    match_newline = re.search(r"\n", line)
                    u_iP_f_ = np.array(line[match_equal.start() + 2:match_newline.start()].split(','), dtype="float64")

                # print "force_name_ =", force_name_
                # print "body_id_ =", body_id_
                # print "Fx_ =", Fx_
                if (force_name_ != None) and (body_id_ != None) and (u_iP_f_ != None) and (Fx_ != None) and (Fy_ != None) and (Mz_ != None):
                    # print "u_iP_f_ =", u_iP_f_
                    #    create force object here
                    force_ = Force(body_id_, force_name=force_name_, Fx=Fx_, Fy=Fy_, Mz=Mz_, u_iP_f=u_iP_f_, parent=parent)
                    # print "object created"
                    forces.append(force_)
                    
                    #    force
                    force_name_ = None
                    body_id_ = None
                    Fx_ = None
                    Fy_ = None
                    Mz_ = None
                    u_iP_f_ = None
        # logging.getLogger("DyS_logger").info("File %s found! Forces created successfully!", filename)
        
    except:
        pass
        # logging.getLogger("DyS_logger").info("File %s not found! No forces created!", filename)
        
    return forces

        
if __name__ == "__main__":
    list = create_list("forces.txt")
    for object in list:
        pprint(vars(object))


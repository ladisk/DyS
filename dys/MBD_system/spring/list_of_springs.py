'''
Created on 18. mar. 2014

@author: lskrinjar (email: skrinjar.luka@gmail.com)
'''
import logging
import os
from pprint import pprint
import re
import time

import numpy as np
from spring import Spring


def create_list(filename, parent=None):
    """
    lists file joints.txt in folder - folder_path and returns the list of joints as objects
    in:
        filename - absolute path to file
    returns:
        springs - list of contact objects
    """
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
            
            for line in file_:
                if line.startswith("SPRINGS FILE-START") or line.startswith("#"):
                    continue
                elif line.startswith("SPRINGS FILE-END"):
    #                     print "Reading file: "+filename+" finished."
                    break
                
                elif line.startswith("spring_name"):
                    match_equal = re.search(r"=", line)
                    match_newline = re.search(r"\n", line)
                    spring_name_ = line[match_equal.start() + 1:match_newline.start()]
                
                elif line.startswith("spring_type"):
                    match_equal = re.search(r"=", line)
                    match_newline = re.search(r"\n", line)
                    spring_type_ = line[match_equal.start() + 1:match_newline.start()]
                
                elif line.startswith("k ="):
                    match_equal = re.search(r"=", line)
                    match_newline = re.search(r"\n", line)
                    k_ = np.float(line[match_equal.start() + 1:match_newline.start()])
            
                elif line.startswith("k_t ="):
                    match_equal = re.search(r"=", line)
                    match_newline = re.search(r"\n", line)
                    k_t_ = np.float(line[match_equal.start() + 1:match_newline.start()])
                    
                elif line.startswith("c ="):
                    match_equal = re.search(r"=", line)
                    match_newline = re.search(r"\n", line)
                    c_ = np.float(line[match_equal.start() + 1:match_newline.start()])
                    
                elif line.startswith("c_t ="):
                    match_equal = re.search(r"=", line)
                    match_newline = re.search(r"\n", line)
                    c_t_ = np.float(line[match_equal.start() + 1:match_newline.start()])
                    
                elif line.startswith("l_0"):
                    match_equal = re.search(r"=", line)
                    match_newline = re.search(r"\n", line)
                    l_0_ = np.float(line[match_equal.start() + 1:match_newline.start()])
    
                elif line.startswith("body_i"):
                    match_equal = re.search(r"=", line)
                    match_newline = re.search(r"\n", line)
                    val_ = line[match_equal.start() + 2:match_newline.start()]
                    try:
                        body_i_ = np.int(val_)
                    except:
                        body_i_ = str(val_)
                    
                elif line.startswith("body_j"):
                    match_equal = re.search(r"=", line)
                    match_newline = re.search(r"\n", line)
                    val_ = line[match_equal.start() + 2:match_newline.start()]
                    try:
                        body_j_ = np.int(val_)
                    except:
                        body_j_ = str(val_)
                    
                    
                elif line.startswith("u_iP"):
                    match_equal = re.search(r"=", line)
                    match_newline = re.search(r"\n", line)
                    u_iP_ = np.array(line[match_equal.start() + 2:match_newline.start()].split(','), dtype="float64")
                    
                elif line.startswith("u_jP"):
                    match_equal = re.search(r"=", line)
                    match_newline = re.search(r"\n", line)
                    u_jP_ = np.array(line[match_equal.start() + 2:match_newline.start()].split(','), dtype="float64")
                    
                if (spring_name_ != None) and (spring_type_ == "translational") and (body_i_ != None) and (body_j_ != None) and (u_iP_ != None) and (u_jP_ != None):
                    #    create spring object here
                    spring_ = Spring(spring_name_, spring_type_, body_i_, body_j_, u_iP_, u_jP_, k=k_, c=c_, l_0=l_0_, parent=parent)
                    springs.append(spring_)
                    
    
                    #    spring type
                    spring_name_ = None
                    spring_type_ = None
                    k_ = 0
                    c_ = 0
                    body_i_ = None
                    body_j_ = None
                    u_iP_ = None
                    u_jP_ = None
                    l_0_ = None
                    
                    
                elif (spring_name_ != None) and (spring_type_ == "torsional") and (body_i_ != None) and (body_j_ != None):
                    #    create spring object here
                    spring_ = Spring(_name=spring_name_, spring_type=spring_type_, k=k_t_, c=c_t_, body_id_i=body_i_, body_id_j=body_j_, parent = parent)
                    springs.append(spring_)
                    
                    #    spring type
                    spring_name_ = None
                    spring_type_ = None
                    k_t_ = 0
                    c_t_ = 0
                    body_i_ = None
                    body_j_ = None

            # logging.getLogger("DyS_logger").info("File %s found! Springs created successfully!", filename)
            file_.close
    else:
        logging.getLogger("DyS_logger").info("File %s not found! No springs created!", filename)
                
    return springs

        
if __name__ == "__main__":
    t1 = time.time()
    spring_list = create_list("springs.txt")
    print time.time() - t1
    for spring in spring_list:
        pprint(vars(spring))
    

    
    

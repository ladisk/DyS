'''
Created on 29. jun. 2013

@author: lskrinjar
'''

import os
import re
from collections import OrderedDict

import numpy as np

from MBD_system.string2array import string2array

def read_data(filename):
    """
    method reads geometry data file of a body
    in:
        filename of the body data file
    out:
        _dict - dictionary of attributes to assgin to body object
    """
    #    predefine empty dictionary
    dict_ = OrderedDict([])

    with open(filename, 'r') as file_:
    #    'r' - read
    #    'w' - write
    #    'a' - append
    #    'r+' - read and write
    
        # #    predefined values
        # name = []
        # body_geometry_filename = []
        #
        # density = 0.
        # volume = 0.
        # mass = 0.
        # J_zz = 0.
        # CM_LCS = np.zeros(3, dtype="float32")
        # R_GCS = np.zeros(3, dtype="float32")
        # theta_GCS = np.zeros(3, dtype="float32")
        # dR_GCS = np.zeros(3, dtype="float32")
        # dtheta_GCS = np.zeros(3, dtype="float32")
        #
        # color_GL = np.ones(3, dtype="float32")
        # transparent_GL = 1
        # display_style = "filled"

        # ln = "name = actuating_lever"
        # print ln[ln.index("=")]

        for line in file_:
            if "#" in line:
                line = line[:line.index("#")].strip()

            #   if line is empty or begins with #-comment, skip it
            if len(line.strip()) == 0 or line.startswith("#"):
                pass
            else:
                #   attribute name
                key = line[0:line.find("=")].strip()

                #   if vasriable name includes units, remove it
                if "[" in key:
                    key = key[0:line.find("[")].strip()


                val = line[line.find("=")+1:].strip()
                

                #   if string has commas, then convert to array
                if "," in val:
                    val = string2array(val)
                else:

                    try:
                        #   try covnert string to integer
                        val = int(val)
                    except:
                        try:
                            #   try convert string to float
                            val = float(val)
                        except:
#                           #    check if string is bool and convert it to bool
                            if val.lower() == "true":
                                val = True
                            elif val.lower() == "false":
                                val = False
                    

                #    add dict item to dict with value
                # print key," = ", val, type(val)
                dict_[key] = val



            # if line.startswith("name="):
            #     line_ = file_.next()
            #     match_newline = re.search(r"\n", line_)
            #     name = line_[:match_newline.start()]
            #
            # elif line.startswith("density"):
            #     density = float(file_.next().strip())
            #
            # elif line.startswith("volume"):
            #     volume = float(file_.next().strip())
            #
            # elif line.startswith("mass"):
            #     mass = float(file_.next().strip())
            #
            # elif line.startswith("J_zz"):
            #     J_zz = float(file_.next().strip())
            #
            # elif line.startswith("CM in CAD LCS"):
            #     CM_CAD_LCS = np.array(file_.next().split(','), dtype="float32")
            #
            # elif line.startswith("CAD LCS in GCS"):
            #     CAD_LCS_GCS = np.array(file_.next().split(','), dtype="float32")
            #
            # elif line.startswith("theta"):
            #     theta_GCS = np.deg2rad(np.array(file_.next().split(','), dtype="float32"))
            #
            # elif line.startswith("dR"):
            #     dR_GCS = np.array(file_.next().split(','), dtype="float32")
            #
            # elif line.startswith("angular"):
            #     dtheta_GCS = np.deg2rad(np.array(file_.next().split(','), dtype="float32"))
            #
            # elif line.startswith("body geometry filename"):
            #     line_ = file_.next()
            #     match_newline = re.search(r"\n", line_)
            #     body_geometry_filename = line_[:match_newline.start()]
            #
            # elif line.startswith("color"):
            #     color_GL = np.array(file_.next().split(','), dtype="float32")
            #
            # elif line.startswith("transparent"):
            #     transparent_GL = float(file_.next().strip())
            #
            # elif line.startswith("display style"):
            #     line_ = file_.next()
            #     match_newline = re.search(r"\n", line_)
            #     display_style = line_.translate(None, '\n')
            #
            # elif line.startswith("material"):
            #     line_ = file_.next()
            #     match_newline = re.search(r"\n", line_)
            #     material = line_[0:match_newline.start()]

    
    #    calculate mass from volume and density
    # if (volume != 0) and (density != 0):
    #     mass = density * volume

    #    close opened file
    file_.closed

    return dict_
    # return name, density, volume, mass, J_zz, CM_CAD_LCS, CAD_LCS_GCS, theta_GCS, dR_GCS, dtheta_GCS, body_geometry_filename, color_GL, transparent_GL, display_style

if __name__ == "__main__":
    # actual_body_name, density, volume, mass, J_zz, CM_CAD_LCS, CAD_LCS_GCS, theta_GCS, dR_GCS, dtheta_GCS, body_geometry_filename, color_GL, transparent_GL, display_style = read_data("data_body_2.txt")
    # print actual_body_name, density, volume, mass, J_zz, CM_CAD_LCS, CAD_LCS_GCS, theta_GCS, dR_GCS, dtheta_GCS, body_geometry_filename, color_GL, transparent_GL, display_style
    dict = read_data("data_body_2.txt")
    for key in dict:
        print key, dict[key], type(dict[key])

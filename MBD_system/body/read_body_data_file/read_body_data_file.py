"""
Created on 29. jun. 2013

@author: lskrinjar
"""

import os
import re
from collections import OrderedDict

import numpy as np

from MBD_system.string2array import string2array

def read_body_data_file(filename):
    """
    method reads geometry data file of a body
    in:
        filename of the body data file
    out:
        _dict - dictionary of attributes to assgin to body object
    """
    #    predefine empty dictionary
    dict_ = OrderedDict([])

    if os.path.isfile(filename):
        #   'r' - read
        #   'w' - write
        #   'a' - append
        #   'r+' - read and write
        with open(filename, 'r') as file_:
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
                            #   try convert string to integer
                            val = int(val)
                        except:
                            try:
                                #   try convert string to float
                                val = float(val)
                            except:
                                #   check if string is bool and convert it to bool
                                if val.lower() == "true":
                                    val = True
                                elif val.lower() == "false":
                                    val = False

                    #    add dict item to dict with value
                    dict_[key] = val

        #    close opened file
        file_.close()
    else:
        print "Body data file not found! Check filename %s"%filename

    return dict_


if __name__ == "__main__":
    dict = read_body_data_file("data_body_2.dat")
    for key in dict:
        print key, dict[key], type(dict[key])

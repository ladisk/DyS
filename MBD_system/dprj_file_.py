'''
Created on 12. avg. 2014

@author: lskrinjar

'''
import os
import ast

def dprj_file_read(filename):
    """
    Method to read file .dprj and save all parameters to dictionary
    Args:
    filename - absolute path to file filename.dprj (dprj - Dynamic PRoJect)
    Raises:
    None
    Returns:
    dict_ - dictionary of all data - properties
    """
    #    predefine empty dictionary
    dict_ = {}
    

    with open(filename, 'r') as file_:
    #    'r' - read
    #    'w' - write
    #    'a' - append
    #    'r+' - read and write
        for line in file_:
            if line.startswith("simulation settings-start"):
                continue

            elif line.startswith("simulation settings-end"):
                #    Close opened file
                file_.closed

            else:
                #    Add a new item to dictionary dict_
                try:
                    value_ = float(line[line.index('=') + 1:].strip())

                except:
                    value_ = str(line[line.index('=') + 1:].strip())

                    if value_.lower() == "true":
                        value_ = True
                    elif value_.lower() == "false":
                        value_ = False

                key_ = line[0:line.index("=")].strip()
                dict_[key_] = value_
        
        return dict_

    
    
def dprj_file_write():
    
    return None


if __name__ == "__main__":
    _dynamic_system = "dynamic_system_0_6_1"
    _dprj_file = "dys_0_6_1.dprj"
    _file = os.path.join(os.path.dirname(os.getcwd()), "dynamic_systems", _dynamic_system, _dprj_file)
    # print os.path.join(os.getcwd(), ..)
    d_ = dprj_file_read(_file)
    print d_

'''
Created on 29. jul. 2015

@author: skrinjar.luka@gmail.com

'''
import numpy as np


def string2array(string):
    """
    
    """
    array = np.array(str(string).strip().split(","), dtype=np.float32)
    return array


if __name__ == "__main__":
    string = " 0.,0."
    print string2array(string)
"""
Created on 29. jul. 2015

@author: skrinjar.luka@gmail.com

"""
import numpy as np


def array2string(array):
    """
    
    """
    string = ""
    for val in array:
        if string == "":
            string = str(val)
        else:
            string = string + ", " + str(val)
    
    return string


if __name__ == "__main__":
    array = np.random.rand(2)
    print array2string(array)
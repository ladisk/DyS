"""
Created on 22. maj 2015

@author: lskrinjar
"""

import math

def convert_size(size_in_bytes):
    """
    Converts size in bytes to kBytes, MBytes,...
    Args:
        size_in_bytes - integer in units of bytes
    """
    size_name = ("B", "KB", "MB", "GB", "TB", "PB", "EB", "ZB", "YB")

    i = int(math.floor(math.log(size_in_bytes, 1024)))
    p = math.pow(1024, i)
 
    s = round(size_in_bytes / p, 2)
     
    if (s > 0):
        return str(s) + size_name[i]
    else:
        return '0B'
    

if __name__ == "__main__":
    print convert_size(10000)

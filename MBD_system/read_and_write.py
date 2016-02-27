"""
Created on 20. feb. 2015

@author: lskrinjar

"""
from pprint import pprint

import dill


def write(obj, filename):
    
    with open(filename, 'wb') as output_:
        dill.dump(obj, output_)


def read(filename):
    with open(filename, 'rb') as input_:
        return dill.load(input_)
    
    
if __name__ == "__main__":
    a = read("ground_.pkl")
    pprint(vars(a))
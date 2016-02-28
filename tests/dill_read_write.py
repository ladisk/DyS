"""
Created on 20. feb. 2015

@author: lskrinjar

"""
import os

import dill


def write(obj, filename):
    os.chdir(os.path.abspath(os.path.join(os.path.dirname(__file__))))
    filename = os.path.join(os.getcwd(), filename)
    
    with open(filename, 'wb') as output_:
        dill.dump(obj, output_)
        print "Object saved to file: ", filename
        

def read(filename):
    with open(filename, 'rb') as input_:
        return dill.load(input_)
    
    

if __name__ == "__main__":
    read("ground_.pkl")
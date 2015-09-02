'''
Created on 20. feb. 2015

@author: lskrinjar

'''
import os
import dill
import pickle
import cPickle

def write(obj, filename):
    
    with open(filename, 'wb') as output_:
        dill.dump(obj, output_)


def read(filename):
    with open(filename, 'rb') as input_:
        return dill.load(input_)
    
    
if __name__ == "__main__":
    read("ground_.pkl")
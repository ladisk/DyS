"""
Created on 20. feb. 2015

@author: lskrinjar

"""
import copy
from pprint import pprint
import pickle
import dill

def read_dill(filename):
    """

    :param filename:
    :return:
    """
    with open(filename, 'rb') as input_:
        return dill.load(input_)


def write_dill(obj, filename):
    """

    :param filename:
    :return:
    """
    with open(filename, 'wb') as output_:

        _obj = copy.copy(obj)
        if hasattr(_obj, "dys"):
            _obj.dys = None
        _obj._parent = None

        _obj._parent = None

        dill.dump(_obj, output_)
        print "Sucessfully saved to file: ", filename


def read_pickle(filename):
    """

    :param filename:
    :return:
    """
    with open(filename, 'rb') as input_:
        return pickle.load(input_)


def write_pickle(obj, filename):
    """

    :param obj:
    :param filename:
    :return:
    """
    with open(filename, 'wb') as output:
        pickle.dump(obj, output, pickle.HIGHEST_PROTOCOL)
        print "Sucessfully saved to file: ", filename


def write_txt(obj, filename):
    """

    :param filename:
    :return:
    """


def read_txt(filename):
    """

    :param filename:
    :return:
    """


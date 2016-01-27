__author__ = 'lskrinjar'

import os

def check_filename(_filename):
    """

    :return:
    """
    #   get filename and extension
    __filename, __extension = os.path.splitext(_filename)
    i = 0
    while os.path.isfile(_filename):
        i += 1
        _filename = __filename + '_%02d'%i + __extension

    return _filename


if __name__ == "__main__":
    _file = "ball.pkl"
    print check_filename(_file)
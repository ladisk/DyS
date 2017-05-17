"""
Created on 11. apr. 2016

@author: skrinjar.luka@gmail.com

"""
import itertools

def list_flatten(lst):
    """
    
    """
    _lst = []
    for item in lst:
        if isinstance(item, list):
            _l = item
        else:
            _l = [item]
        
        _lst.append(_l)

    return list(itertools.chain(*_lst))

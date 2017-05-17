'''
Created on 21. maj 2015

@author: lskrinjar

'''
import numpy as np
import copy
import time

def gaussian_elimination(A, b, type=None):
    """
    Function returns a matrix based on input matrix on which a gaussian
    elimination has been performed
    Args:
        A - matrix
        b - vector (optional)
        type - type of pivoting
             - None
             - PP
             - SPP
             - FP
    """
    #    augmented matrix - matrix A augmented with vector b
    if b is not None:
        A = np.column_stack((A, b))
    
    print A
    i_rows, j_cols = np.shape(A)


    col_indx = np.arange(0, j_cols)
    row_indx = np.arange(0, i_rows)


    for i in range(0, i_rows - 1):
        #    i-th row of matrix A
        i_row = A[i, :]

        #    sub column of diagonal element
        abs_sub_col = abs(A[i + 1:, i])
        
        #    max element of sub column
        a_max_sub_col = max(abs_sub_col)
        pos_max_sub_col = abs_sub_col.argmax()

        if a_max_sub_col == 0:
            print A
            print "i-th row =", i
            print "row =", i_row
            print "sub col =", A[i + 1:, i]
            raise ValueError("Matrix is singular.")

        
        #    row multiplier
        m = A[i, i] / a_max_sub_col
        
        #    copy rows to replace the in the next two steps
        _i_row_indx = row_indx[i]
        _max_row_indx = i + pos_max_sub_col + 1
        
        _A_i_row = copy.copy(i_row)
        _A_max_row = A[i + pos_max_sub_col + 1, :]
        
        
        #    replace rows
        row_indx[i] = _max_row_indx
        row_indx[i + pos_max_sub_col + 1] = _i_row_indx
        
        A[i, :] = _A_max_row
        A[i + pos_max_sub_col + 1, :] = _A_i_row

        for i_sub_row in range(i + 1, i_rows):
            #    row k - row i
            A[i_sub_row, :] = A[i_sub_row, :] - m * i_row
    print "GE-start"
    print A 
    print row_indx
    print "GE-end"
    return None

if __name__ == "__main__":
    A = np.array([[0.003, 59.14],
                  [5.291, -6.13]])
    b = np.array([59.17, 46.78])
    sol = gaussian_elimination(A, b)
    print "sol ="
    print sol

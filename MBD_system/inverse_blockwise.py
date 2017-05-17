'''
Created on 22. apr. 2015

@author: luka.skrinjar
'''
import time
import numpy as np
np.set_printoptions(precision=4, threshold=None, edgeitems=100, linewidth=1000, suppress=False, nanstr=None, infstr=None)


def inverse_blockwise(A_inv, B, C, D, Adim):
    """
    Link: http://en.wikipedia.org/wiki/Invertible_matrix
    m = [A B
         C D]
         
    Args:
        matrix - 
        dim_sub_mat_A - 
        M - mass matrix
    """
#    A = matrix[0:dim_sub_mat_A, 0:dim_sub_mat_A]
#    B = matrix[0:dim_sub_mat_A, dim_sub_mat_A:]
#    C = matrix[dim_sub_mat_A:, 0:dim_sub_mat_A]
#    D = matrix[dim_sub_mat_A:, dim_sub_mat_A:]

    # print "M_inv =", A_inv
    # print "B =", B
    # print "C =", C
    # print "D =", D
    # time.sleep(100)
    #    inverse of submatrix A
#     print "C =", C.shape
#     print "A_inv =", A_inv.shape
#     print "B =", B.shape
#     print "reduce(np.dot, [C, A_inv, B]) =", reduce(np.dot, [C, A_inv, B])
    #    inverse of sub matrix
    D_C_A_inv_B__inv = np.linalg.inv(D - reduce(np.dot, [C, A_inv, B]))

    #    calculating sub matrices of inverse matrix
    _matrix_11 = +A_inv + reduce(np.dot, [A_inv, B, D_C_A_inv_B__inv, C, A_inv])
    _matrix_12 = reduce(np.dot, [-A_inv, B, D_C_A_inv_B__inv])
    _matrix_22 = D_C_A_inv_B__inv

    #    is submatrix D is zero matrix (when Cs=0), the submatrix _matrix21 equals the transpose of _matrix12*(-1)
    if np.count_nonzero(D) == 0:
        _matrix_21 = -1 * _matrix_12.T
    else:
        _matrix_21 = reduce(np.dot, [-D_C_A_inv_B__inv, C, -A_inv])

    #    insert calculated inverse submatrices in predefined zero matrix
    _matrix = np.empty([np.shape(A_inv)[1] + np.shape(B)[1], np.shape(A_inv)[0] + np.shape(C)[0]])
    
    _matrix[0:Adim, 0:Adim] = _matrix_11
    _matrix[0:Adim, Adim:] = _matrix_12
    _matrix[Adim:, 0:Adim] = _matrix_21
    _matrix[Adim:, Adim:] = _matrix_22

    return _matrix


if __name__ == "__main__":
    m = np.array([[1, 0, 0, 0, 5, 8],
                  [0, 1, 0, 0, 0, 3],
                  [0, 0, 0, 0, 1, 2],
                  [0, 0, 0, 0, 5, 1],
                  [5, 0, 1, 5, 0, 0],
                  [8, 3, 2, 1, 0, 0]])

#    m = np.array([[1, 0, 2, 1],
#                  [0, 1, 0, 2],
#                  [2, 0, 0, 0],
#                  [1, 2, 0, 0]])

    t0 = time.time()
    a = inverse_blockwise(m, 2)
#    a = np.linalg.inv(m)
    print "time =", time.time() - t0
    print a
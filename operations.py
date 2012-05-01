#!/usr/bin/env python
import math

"""
Module to accompany genstruct.  This will contain the operations needed
to reproduce numpy operations.  The idea is to eventually phase out
numpy from genstruct.
"""

zeros3 = [[0.,0.,0.],[0.,0.,0.],[0.,0.,0.]]
zeros1 = [0.,0.,0.]
zeros6 = [0.,0.,0.,0.,0.,0.]
identity3 = [[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]]
def vect_length(vect1, vect2):
    """
    Gives the distance between two vectors
    """
    dvect = vect_sub(vect1, vect2)
    dprod = dot(dvect, dvect)
    return math.sqrt(dprod)

def vect_add(vect1, vect2):
    """
    Adds elements of each vector together.
    """
    return [i+j for i, j in zip(vect1, vect2)]

def vect_sub(vect1, vect2):
    """
    Subtracts elements of each vector together.
    """
    return [i-j for i, j in zip(vect1, vect2)]

def matrx_sub(mat1, mat2):
    """
    Subtracts one matrix from another.
    """
    return [vect_sub(i,j) for i, j in zip(mat1, mat2)]

def scalar_mult(scalar, vector):
    """
    Returns a scalar multiple of the vector.
    """
    return [i*scalar for i in vector]

def scalar_div(scalar, vector):
    """
    Returns a scalar quotient of the vector.
    """
    return [i/scalar for i in vector]

def matrx_mult(vector, matrix):
    """
    Returns a vector resulting from multiplication of a vector
    with a matrix.
    """
    #TODO(pboyd): add tests to see if the matrix and vector
    # can be multiplied.
    transmat = transpose(matrix)
    return [dot(vector, i) for i in transmat]

def transpose(matrix):
    """Returns a transpose of a matrix."""
    #TODO(pboyd): make sure the matrix is square
    return [list(i) for i in zip(*matrix)]

def dot(vector1, vector2):
    """
    returns the dot product of two vectors.
    """
    return sum([x*y for x,y in zip(vector1, vector2)])

def inner(vector1, vector2):
    """
    Returns the inner product of two vectors.
    """
    return [x*y for x,y in zip(vector1, vector2)]

def cross(vector1, vector2):
    """
    Returns the cross produc of two vectors.
    3-D only for now!!
    """
    v0 = (vector1[1] * vector2[2]) - (vector1[2] * vector2[1])
    v1 = -1.*((vector1[0] * vector2[2]) - (vector1[2] * vector2[0]))
    v2 = (vector1[0] * vector2[1]) - (vector1[1] * vector2[0])
    return [v0, v1, v2]


#!/usr/bin/env python
import math
from scipy.spatial import distance

"""
Module to accompany genstruct.  This will contain the operations needed
to reproduce numpy operations.  The idea is to eventually phase out
numpy from genstruct.
"""

# Constants

RAD2DEG = 180./math.pi
DEG2RAD = math.pi/180.
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

def bond_tolerance(atom1, atom2):
    """
    returns a tolerance value for atom - atom distances
    case specific
    """
    # This function will need some elaboration
    # TODO(pboyd): these conditions are not robust and
    # should be changed when the code becomes bigger
    # G is the metal - metal bond for linked chains
    if("G" in (atom1, atom2) and 
        (("Y" not in (atom1, atom2))and("Z" not in (atom1, atom2)))):
        return 1.35
    elif("X" in (atom1, atom2) and 
        (("Y" not in (atom1, atom2))and("Z" not in (atom1, atom2)))):
        return 0.97
    elif("G" in (atom1, atom2) or "X" in (atom1, atom2))and \
        (("Y" in (atom1, atom2))or("Z" in (atom1, atom2))):
        return 1.1
    elif("Y" in (atom1, atom2))and("G" not in (atom1, atom2) 
            or "X" not in (atom1,atom2)):
        return 0.
    elif("Z" in (atom1, atom2))and("G" not in (atom1, atom2) 
            or "X" not in (atom1,atom2)):
        return 0.
    elif("Br" in (atom1, atom2))and("J" in (atom1, atom2)):
        return 2.0
    elif("I" in (atom1, atom2))and("J" in (atom1, atom2)):
        return 2.5
    elif("Cl" in (atom1, atom2))and("J" in (atom1, atom2)):
        return 1.8
    elif("H" in (atom1, atom2))and("C" in (atom1, atom2)):
        return 1.15
    elif("P" in (atom1, atom2))and("C" in (atom1, atom2)):
        return 1.9
    elif(set((atom1, atom2)) == set("C")):
        return 1.6
    elif("Cu" in (atom1, atom2))and("N" in (atom1,  atom2)):
        return 2.0
    elif("Zn" in (atom1, atom2))and("N" in (atom1, atom2)):
        return 2.0
    elif("In" in (atom1, atom2))and("O" in (atom1, atom2)):
        return 2.0
    elif("V" in (atom1, atom2))and("O" in (atom1, atom2)):
        return 2.0
    elif("Ba" in (atom1, atom2))and("O" in (atom1, atom2)):
        return 2.6
    elif("O" in (atom1, atom2)):
        return 0.
    elif(set((atom1, atom2)) == set("Ba")):
        return 4.6
    elif(set((atom1, atom2)) == set("H")):
        return 0.
    elif("I" in (atom1, atom2))and("C" in (atom1, atom2)):
        return 2.2
    elif("H" in (atom1, atom2)):
        return 0.
    else:
        return 0. 

def length(coord1, coord2=None):
    """ 
    Returns the length between two vectors.
    If only one vector is specified, it returns
    the length of that vector from the origin
    """
    if coord2 is not None:
        coord2 = coord2
    else:
        coord2 = zeros1[:]
    dist = distance.cdist([coord1], [coord2], 'euclidean')[0][0]
    return dist 

def calc_angle(vect1, vect2):
    """ determines angle between vector1 and vector2"""
    dot12 = dot(vect1, vect2) 
    dist11 = length(vect1)
    dist22 = length(vect2)
    # clamps the angle coefficient to a min or max of -1, 1 so no 
    # error is returned when calculating the acos.
    angle_coefficient =  min(max(dot12/(dist11*dist22), -1.0),1.0)
    return math.acos(angle_coefficient)


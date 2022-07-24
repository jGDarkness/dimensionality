# Import the necessary libraries.
import os
import scipy.special
from math import comb

class shape():
    def points(dimensions):
        points = 2**dimensions
        points = "{:,}".format(points)
        return points
    
    def edges(dimensions):
        edges = dimensions * (2**(dimensions - 1))
        edges = "{:,}".format(edges)
        return edges
    
    def square(dimensions):
        dimensions = dimensions - 1 
            # For some reason, this function returns the results for one higher dimension than requested. Need to research why this is.
        lines = (dimensions) * (dimensions + 1) * (2**(dimensions - 2))
        lines = "{:,}".format(lines)
        return lines

    def cube(dimensions):
        cube = scipy.special.binom(dimensions, 3) * 2**(dimensions - 3)
        cube = "{:,}".format(cube)
        return cube

    def tesseract(dimensions):
        tesseract = 2**(dimensions - 4) * comb(dimensions, 4)
        tesseract = "{:,}".format(tesseract)
        return tesseract
    
    def penteract(dimensions):
        penteract = 2**(dimensions-5) * scipy.special.binom(dimensions, 5)
        penteract = "{:,}".format(penteract)
        return penteract

    def hexeract(dimensions):
        dimensions = dimensions - 6
            # For some reason, this function returns the results for 6 higher dimensions than requested. Need to research why this is.
        hexeract = (2**dimensions) * comb(dimensions+6,6)
        hexeract = "{:,}".format(hexeract)
        return hexeract

    def hepteract(dimensions):
        hepteract = 2**(dimensions-7)*scipy.special.binom(dimensions, 7)
        hepteract = "{:,}".format(hepteract)
        return hepteract

    def octeract(dimensions):
        dimensions = dimensions - 8
            # For some reason, this function returns the results for 8 higher dimensions than requested. Need to research why this is.
        octeract = scipy.special.binom(dimensions+8,8) * (2**dimensions)
        octeract = "{:,}".format(octeract)
        return octeract

    def enneract(dimensions):
        dimensions = dimensions - 9
            # For some reason, this function returns the results for 9 higher dimensions than requested. Need to research why this is.
        enneract = scipy.special.binom(dimensions+9,9) * (2**dimensions)
        enneract = "{:,}".format(enneract)
        return enneract

    def dekeract(dimensions):
        dimensions = dimensions - 10
            # For some reason, this function returns the results for 9 higher dimensions than requested. Need to research why this is.
        dekeract = scipy.special.binom(dimensions+9,9) * (2**dimensions)
        dekeract = "{:,}".format(dekeract)
        return dekeract

    def hypercube(dimensions):
        originalDimensions = dimensions
        dimensions = dimensions - dimensions
            # For some reason, this function returns the results for 9 higher dimensions than requested. Need to research why this is.
        hypercube = scipy.special.binom(dimensions+originalDimensions-1,originalDimensions-1) * (2**dimensions)
        hypercube = "{:,}".format(hypercube)
        return hypercube

# number E_sub(m,n) of m-dimensional faces of a n-dimensional hypercube
#                        0face  1face  2face  3face  4face 5face 6face 7face 8face 9face 10face
pointArray =            (    1,     0,     0,     0,     0,    0,    0,    0,    0,    0,     0) # 0D
lineArray =             (    2,     1,     0,     0,     0,    0,    0,    0,    0,    0,     0) # 1D
squareArray =           (    4,     4,     1,     0,     0,    0,    0,    0,    0,    0,     0) # 2D
cubeArray =             (    8,    12,     6,     1,     0,    0,    0,    0,    0,    0,     0) # 3D
tesseractArray =        (   16,    32,    24,     8,     1,    0,    0,    0,    0,    0,     0) # 4D
penteractArray =        (   32,    80,    80,    40,    10,    1,    0,    0,    0,    0,     0) # 5D
hexeractArray =         (   64,   192,   240,   160,    60,   12,    1,    0,    0,    0,     0) # 6D
hepteractArray =        (  128,   448,   672,   560,   280,   84,   14,    1,    0,    0,     0) # 7D
octeractArray =         (  256,  1024,  1792,  1792,  1120,  448,  112,   16,    1,    0,     0) # 8D
enneractArray =         (  512,  2304,  4608,  5376,  4032, 2016,  672,  144,   18,    1,     0) # 9D
dekeractArray =         ( 1024,  5120, 11520, 15360, 13440, 8064, 3360,  960,  180,   20,     1) # 10D

# The number of "0-face"  Points (Vertex)                in each Array[1]  is 2^n                                 
# The number of "1-face"  Lines Segments (Edge)          in each Array[2]  is n*2^(n-1)                           
# The number of "2-face"  Squares (Face)                 in each Array[3]  is where n = n-1, let n*(n+1)*2^(n-2)    offset 1 zeroes
# The number of "3-face"  Cubes (Cell)                   in each Array[4]  is [Binomial[n, 3]*2^(n-3)             
# The number of "4-face"  Tesseracts (hypercube)         in each Array[5]  is 2^(n-4)*C(n,4)                      
# The number of "5-face"  Penteracts (5D hypercubes)     in each Array[6]  is 2^(n-5)*binomial(n,5)               
# The number of "6-face"  Hexeracts (n+6D hypercubes)    in each Array[7]  is where n = n-5, let 2^n*C(n+6,6)       offset 6 zeroes
# The number of "7-face"  Hepteracts (n+7D hypercubes)   in each Array[8]  is 2^(n-7)*binomial(n,7)               
# The number of "8-face"  Octeracts (n+8D hypercubes)    in each Array[9]  is binomial(n+8,8) * 2^n                 offset 8 zeroes
# The number of "9-face"  Enneracts (n+9D hypercubes)    in each Array[10] is binomial(n+9,9) * 2^n                 offset 9 zeroes
# The number of "10-face" Dekeract (n+10D hypercubes)    in each Array[11] is binomial(0+n-1,n-1) * 2^0             offset 10 or more zeroes based on dimensionality

os.system('CLS')

dims = int(input("How many dimensions?\n"))
print("\nA shape with ", dims, " dimensions has:\n- ", 
    shape.points(dims), " Vertex(eces) or Point(s).\n- ",
    shape.edges(dims), " Edge(s) or Line Segment(s).\n- ",
    shape.square(dims), " Face(s) or Square(s).\n- ",
    shape.cube(dims), " Cell(s) or Cube(s).\n- ",
    shape.tesseract(dims), " Tesseract(s) or Hypercube(s).\n- ",
    shape.penteract(dims), " Penteract(s) or 5D Hypercube(s).\n- ",
    shape.hexeract(dims), " Hexeract(s) or 6D Hypercube(s).\n- ",
    shape.hepteract(dims), " Hepteract(s) or 7D-Hypercube(s).\n- ",
    shape.octeract(dims), " Octeract(s) or 8D Hypercube(s).\n- ",
    shape.enneract(dims), " Enneract(s) or 9D Hypercube(s).\n- ",
    shape.dekeract(dims), " Dekeract(s) or 9D Hypercube(s).\n- ",    
)
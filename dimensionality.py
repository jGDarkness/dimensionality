# Libraries
import os
from xmlrpc.client import FastParser
import pyvista                                   # Used to allow interaction with the console and file system.
import scipy.special                        # Used to bring in 2-factor binomial transformation.
from math import comb                       # Used to bring in the mathematical Combination (C) operation.
import numpy as np                          # Used to bring in matrix notation and operations.
import matplotlib.pyplot as plt             # Used to bring in visualization libraries.                    
from mpl_toolkits.mplot3d import Axes3D     
import pyvista as pv
from pyvista import examples


class shape():

    def points(dimensions):
        points = 2**dimensions
        # points = "{:,}".format(points)
            # this causes "invalid literal for int() with base 10: '1,024' (meaning that the comma in the format is breaking the code.
            # If the value of points is to be used elsewhere, this format function should be used at that point in time via another function, 
            # and not part of the class.
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

class matrices():
    def projectionMatrix(dimensions):
        # Should result in a projection matrix one row less than the number of dimensions requiested. This is what helps the matrix multiplication
        # reduce the number of coordinates.
        rows = int(0)
        cols = int(0)

        projectionMatrix = np.ndarray(shape=(dimensions-1, dimensions), dtype=np.int32)

        for rows in range(dimensions-1):
            for cols in range(dimensions):
                if cols == rows:
                    projectionMatrix[rows, cols] = 1
                else:
                    projectionMatrix[rows, cols] = 0

        return projectionMatrix

    def verticesMatrix(dimensions, magnitude):
        if dimensions == 0:
            verticesMatrix = np.array([0, 0, 0])
                # Ensures that the return value is a 3D coordinate for plotting. 
            rows = 1
            return verticesMatrix, rows
        else:
            verticesMatrix = magnitude * (2*((np.arange(2**dimensions)[:,None] & (1 << np.arange(dimensions))) > 0) - 1)
            rows, cols = verticesMatrix.shape
                # Handles 3D and all higher dimensional requests.

            if dimensions == 1:
                remainingAxes = np.array([0,0])
                verticesMatrix = np.insert(verticesMatrix,1,remainingAxes,axis=1)
                verticesMatrix = np.insert(verticesMatrix,2,remainingAxes,axis=1)
                # If the user has requested 1D, Ensures that the return value is a 3D coordinate for plotting. 
            elif dimensions ==2:
                remainingAxes = np.array([0,0,0,0])
                verticesMatrix = np.insert(verticesMatrix,2,remainingAxes,axis=1)          
                # If the user has requested 2D, Ensures that the return value is a 3D coordinate for plotting. 

            totalRequestedVertices = int(shape.points(dimensions))

            if rows != totalRequestedVertices:
                print("array doesn't match request")
            else:
                return verticesMatrix, rows
       
def showDimensionality():
    dims = int(input("How many dimensions?\n"))
        # Prompt user for number of dimensions to calculate attributes for.
    magnitude = int(input("What magnitude do you want?\n"))
        #magnitude = int(input("What magnitude do you wish to apply to the shapes?\n"))

    if dims > 2:
        print ("Matrix multiplication has not yet been implemented. You are limited to 0D (points), 1D (lines), and 2D (faces) shapes.")
    else:
        myVertices = matrices.verticesMatrix(dims, magnitude)

        pl = pyvista.Plotter()
        pl.set_scale(xscale=10, yscale=10, zscale=10)
        pl.add_mesh(pyvista.PolyData(myVertices[0], force_float=False))
        pl.view_isometric()
        pl.add_bounding_box()
        pl.show()
    

os.system('CLS')
showDimensionality()
    # if dimensions higher than 3 are requested, conduct matrix multplication on each row of the vertices matrix to generate a 3D vector, and then display.
### LIBRARIES #######################################################################################################################################################
#
import os
from attr import NOTHING
from matplotlib import projections
from pyparsing import null_debug_action                                           
    # Console and file system interation.
import scipy.special               
    # 2-factor binomial transformation. 
from math import comb                   
    # Mathematical "combination" nCr opartion.  
import numpy as np                          
    # Mostly for Matrix/Array notation and operations.
import pyvista as pv                       
from pyvista import examples              
    # Data visulizations
import warnings               
    # Capture and handle warnings.
import tkinter                    
from vtkmodules.tk.vtkTkRenderWindowInteractor import vtkTkRenderWindowInteractor
from zmq import NULL
    # Support for GUI interface using Tk          
#
### END LIBRARIES ###################################################################################################################################################

warnings.filterwarnings('ignore', '.*force_float.*')
    # Suppress 'force_float=True' warning messages from pyvista.

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
    
    def rotationMatrix(rotationAxis, theta):
        rotationXTheta = np.array(
            [1,             0,                0],
            [0, np.cos(theta), -(np.sin(theta))],
            [0, np.sin(theta),    np.cos(theta)]
        )

        rotationYTheta = np.array(
            [   np.cos(theta), 0,  np.sin(theta)],
            [               0, 1,              0],
            [-(np.sin(theta)), 0,  np.cos(theta)]
        )

        rotationZTheta = np.array(
            [np.cos(theta), -(np.sin(theta)), 0],
            [np.sin(theta),    np.cos(theta), 0],
            [            0,                0, 1]
        )

class flatten():
    def myDimension(vertices, dims, magnitude, originalVertices):
        
        myProjection = matrices.projectionMatrix(dims)
        myProjectionRows, myProjectionCols = np.shape(myProjection)

        myVerticesRows, myVerticesColsAsRows = originalVertices[0].shape
        
        if myProjectionCols != myVerticesColsAsRows:
            print ("There was a fatal error in the matrix at dimension:", dims)
        else:
            os.system("CLS")
            row = 0
            print("\nFlattening dimension:", dims)
            myFlattenedVertices = np.ndarray(shape=(myVerticesRows, myProjectionCols-1))
            while row <= vertices - 1 :
                myFlattenedVertices[row] = np.matmul(myProjection, originalVertices[0][row])
                row += 1
            print ("Flattened Vertices Returned:\n", myFlattenedVertices)
            return myFlattenedVertices

def showDimensionality():
    dims = 5
    magnitude = int(1)
        #magnitude = int(input("What magnitude do you wish to apply to the shapes?\n"))
    myVertices = matrices.verticesMatrix(dims, magnitude)
    myVerticesRows, myVerticesColsAsRows = myVertices[0].shape

    myNewVertices = NULL

    if dims > 3:
        dimensions = dims
        
        while dimensions > dims -3:
            myNewVertices = flatten.myDimension(myVerticesRows, dimensions, magnitude, myVertices)
            myVertices[0] = myNewVertices
            dimensions -= 1

        """
        if dimensions == 11:
            myNewVertices = flatten.myDimension(myVerticesRows, dimensions, magnitude, myVertices)
            dimensions = dimensions - 1

        if dimensions == 10:
            myNewVertices = flatten.myDimension(myVerticesRows, dimensions, magnitude, myVertices)
            dimensions = dimensions - 1

        if dimensions == 9:
            myNewVertices = flatten.myDimension(myVerticesRows, dimensions, magnitude, myVertices)
            dimensions = dimensions - 1

        if dimensions == 8:
            myNewVertices = flatten.myDimension(myVerticesRows, dimensions, magnitude, myVertices)
            dimensions = dimensions - 1

        if dimensions == 7:
            myNewVertices = flatten.myDimension(myVerticesRows, dimensions, magnitude, myVertices)
            dimensions = dimensions - 1

        if dimensions == 6:
            myNewVertices = flatten.myDimension(myVerticesRows, dimensions, magnitude, myVertices)
            dimensions = dimensions - 1

        if dimensions == 5:
            myNewVertices = flatten.myDimension(myVerticesRows, dimensions, magnitude, myVertices)
            dimensions = dimensions - 1

        if dimensions == 4:
            myNewVertices = flatten.myDimension(myVerticesRows, dimensions, magnitude, myVertices)

        myVertices = myNewVertices
        if myVertices == None:
            print ("No dimensions were flattened all the way to 3D for projection.")
            quit()
        
        print (myVertices)
        """
    else:
        myVertices = matrices.verticesMatrix(dims, magnitude)

    numberOfPoints = shape.points(dims)
    
    myLabels = ()
    for rows in range(numberOfPoints):
        myLabels = np.append(myLabels, rows+1)
        # Create a list of sequentially numbered labels corresponding to the number of vertices in the final result.

    pl = pv.Plotter()
        
    pl.add_camera_orientation_widget(animate=True, n_frames=180)
    pl.add_axes (line_width=3, labels_off=False)
    pl.add_floor(face='-z', i_resolution=1024, j_resolution=1024, color='black', line_width=3, edge_color='white', opacity=0.2)

    pl.add_mesh(pv.PolyData(myVertices[0]))
    #pl.add_bounding_box(line_width=1, color='white')
    
    """ I need to figure out how to do the unit distance testing between coordinates.
    pl.add_point_labels(
        myVertices[0], 
        myLabels, 
        font_size=12, 
        point_color='red', 
        point_size=15, 
        render_points_as_spheres=True, 
        always_visible=True,
        fill_shape=False)
            # Need to find a way to fix the force_float=False warning. 
        
    pl.add_lines(myVertices[0], color='white', width=2)
    """

    pl.show()

showDimensionality()
    # if dimensions higher than 3 are requested, conduct matrix multplication on each row of the vertices matrix to generate a 3D vector, and then display.
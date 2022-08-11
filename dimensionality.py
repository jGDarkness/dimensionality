### LIBRARIES #######################################################################################################################################################
#
import os
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
from vtkmodules.tk.vtkTkRenderWindowInteractor import vtkTkRenderWindowInteractor
    # Support for GUI interface using Tk          
#
### END LIBRARIES ###################################################################################################################################################

warnings.filterwarnings('ignore', '.*force_float.*')
    # Suppress 'force_float=True' warning messages from pyvista.
warnings.filterwarnings('ignore', '.*auto_close*')
    # Supress 'auto-close' warning from pyvista.

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
            print ("There was a fatal error in the matrix at dimension", dims)
        else:
            os.system("CLS")
            row = 0
            myFlattenedVertices = np.ndarray(shape=(myVerticesRows, myProjectionCols-1))
            while row <= vertices - 1 :
                myFlattenedVertices[row] = np.matmul(myProjection, originalVertices[0][row])
                row += 1
            return myFlattenedVertices, myVerticesRows

def showDimensionality():
    dims = 0
    magnitude = int(1)
        #magnitude = int(input("What magnitude do you wish to apply to the shapes?\n"))
    
    if dims > 3:
        myVertices = matrices.verticesMatrix(dims, magnitude)
        myVerticesRows, myVerticesColsAsRows = myVertices[0].shape
        myNewVertices = np.ndarray(shape=(myVerticesRows, myVerticesColsAsRows), dtype=np.int32)

        pass
        # I can't quite figure out how to recurse through the array flattening process and get an actual result.
        #dimensions = dims
        
        #while dimensions > 3:
        #    myNewVertices = flatten.myDimension(myVerticesRows, dimensions, magnitude, myVertices)
        #    dimensions -= 1
        ###########################################################################################################
    else:
        myNewVertices = matrices.verticesMatrix(dims, magnitude)

    numberOfPoints = shape.points(dims)
    
    myLabels = ()
    for rows in range(numberOfPoints):
        myLabels = np.append(myLabels, rows + 1)
        # Create a list of sequentially numbered labels corresponding to the number of vertices in the final result.

    match dims:
        case 0:
            myLegend = "0D: A Point"
        case 1:
            myLegend = "1D: A Line"
        case 2:
            myLegend = "2D: A Plane"
        case 3:
            myLegend = "3D: A Cube"
        case 4:
            myLegend = "4D: A Tesseract, Flattened to 3D"
        case 5:
            myLegend = "5D: A Penteract, Flattened to 3D"
        case 7:
            myLegend = "6D: A Hexeract, Flattened to 3D"
        case 8:
            myLegend = "7D: An Hepteract, Flattened to 3D"
        case 9:
            myLegend = "8D: An Octeract, Flattened to 3D"
        case 10:
            myLegend = "9D: An Enneract, Flattened to 3D"
        case 11:
            myLegend = "10D: A Tesseract of n=10 Dimensions"

    pl = pv.Plotter()
    
    viewup = [0,0,1]    
    orbit = pl.generate_orbital_path(factor=4.0, n_points=360, shift=3.0, viewup=viewup)
    
    pl.add_axes (line_width=3, labels_off=False)
    pl.add_mesh(pv.PolyData(myNewVertices[0]), smooth_shading=True, label=myLegend, color='white')
    pl.add_point_labels(
        myNewVertices[0], 
        myLabels, 
        font_size=12, 
        point_color='red', 
        point_size=10, 
        render_points_as_spheres=True, 
        always_visible=True,
        fill_shape=False)
    pl.add_legend(face=None)
    pl.add_title('Dimensionality', font='arial', color='white', font_size=16)    

    # ADD CONNECTING LINES

    if dims == 0:
        pass
            # Don't draw any connecting lines.

    elif dims == 1:     # the second coordinate has the same value added to i as 'dim'
        i = 0
        lines = np.array((myNewVertices[0][i], myNewVertices[0][i+1]))
        pl.add_lines(lines)
            # Draw one connecting line between points.

    elif dims == 2:
        i = 0
        while i < 4: # (dim * 2)
            lines = np.array((myNewVertices[0][i], myNewVertices[0][(i + 2) % 2])) # ((i + dim) % dim)
            pl.add_lines(lines)
            i += 1 # (dim / 2)
        i = 0
        while i < 4: # (dim * 2)
            lines2 = np.array((myNewVertices[0][i], myNewVertices[0][(i + 1) % 4])) # (i + (dim / 2)) % (dim / 2))
            pl.add_lines(lines2)
            i += 2 # (dim)
            # Draw four lines connecting the plane.
    
    elif dims == 3:
        i = 0
        while i < 4: # (dim - 1)
            lines = np.array((myNewVertices[0][i], myNewVertices[0][(i+2)%4])) # (i + (dim - 1) % (dim + 1))
            pl.add_lines(lines)
            lines2 = np.array((myNewVertices[0][i+4], myNewVertices[0][((i+2)%4)+4])) # ((i + (dim - 1)) % (dim + 1) + (dim + 1)) 
            pl.add_lines(lines2)
            lines3 = np.array((myNewVertices[0][i], myNewVertices[0][(i+4)])) # (i + (dim + 1))
            pl.add_lines(lines3)
            i += 1 # (dim / 3)
        i = 0
        while i < 4: # (dim - 1)
            lines4 = np.array((myNewVertices[0][i], myNewVertices[0][(i+1)%4])) # (i + (dim - 3) % (dim + 1))
            pl.add_lines(lines4)
            lines5 = np.array((myNewVertices[0][i+4], myNewVertices[0][((i+1)%4)+4])) # ((i + (dims - 3) % (dim + 1)) + (dim + 1))
            pl.add_lines(lines5)
            i += 2 # (dim - 1)
            # Draw lines connecting two planes to each other.
    
    elif dims == 4:
        i = 0
        while i < 4: # (dim)
            lines = np.array((myNewVertices[0][i], myNewVertices[0][(i + 4) % 4])) # ((i + (dim)) % dim)
            pl.add_lines(lines)
            i += 1 #(dim / dim)

    pl.open_movie('orbit.mp4')
    pl.orbit_on_path(orbit, write_frames=True, step=0.0027, viewup=viewup)
    pl.show(auto_close=False)

showDimensionality()
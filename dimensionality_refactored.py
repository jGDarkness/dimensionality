### IMPORTS ###
import scipy.special
from math import comb
import numpy as np
import pyvista as pv
from pyvista import examples
import warnings
from vtkmodules.tk.vtkTkRenderWindowInteractor import vtkTkRenderWindowInteractor
import logging
import threading
from pyvistaqt import BackgroundPlotter
from pyvistaqt import QtInteractor, MainWindow
import sys
import os
from qtpy import QtWidgets
import traceback

### LOGGING ###
logging.basicConfig(filename='dimensionality.log', filemode='w', format='%(asctime)s - %(name)s - %(levelname)s - %(message)s', level=logging.DEBUG)
logging.info("Logging has been configured.")

def exc_handler(exctype, value, tb):
    logging.exception("".join(traceback.format_exception(exctype, value, tb)))
sys.excepthook = exc_handler
logging.info("Traceback logging has been enabled.")

### WARNINGS HANDLERS ###
warnings.filterwarnings('ignore', '.*force_float.*')
warnings.filterwarnings('ignore', '.*auto_close*')
logging.info("Warning filters enabled: '.*force_float.*' and '.*auto_close*'.")

os.environ["QT_API"] = "pyqt5"
logging.info("Setting 'QT_API' to 'pyqt5'.")

class shape():    
    def points(dimensions):
        logging.info("Class.def 'shape.points' invoked.")
        points = 2**dimensions
        return points
    
    def edges(dimensions):
        logging.info("Class.def 'shape.edges' invoked.")
        edges = dimensions * (2**(dimensions - 1))
        edges = "{:,}".format(edges)
        return edges
    
    def square(dimensions):
        logging.info("Class.def 'shape.square' invoked.")
        dimensions = dimensions - 1 
            # For some reason, this function returns the results for one higher dimension than requested. Need to research why this is.
        lines = (dimensions) * (dimensions + 1) * (2**(dimensions - 2))
        lines = "{:,}".format(lines)
        return lines

    def cube(dimensions):
        logging.info("Class.def 'shape.cube' invoked.")
        cube = scipy.special.binom(dimensions, 3) * 2**(dimensions - 3)
        cube = "{:,}".format(cube)
        return cube

    def tesseract(dimensions):
        logging.info("Class.def 'shape.tesseract' invoked.")
        tesseract = 2**(dimensions - 4) * comb(dimensions, 4)
        tesseract = "{:,}".format(tesseract)
        return tesseract
    
    def penteract(dimensions):
        logging.info("Class.def 'shape.penteract' invoked.")
        penteract = 2**(dimensions-5) * scipy.special.binom(dimensions, 5)
        penteract = "{:,}".format(penteract)
        return penteract

    def hexeract(dimensions):
        logging.info("Class.def 'shape.hexeract' invoked.")
        dimensions = dimensions - 6
            # For some reason, this function returns the results for 6 higher dimensions than requested. Need to research why this is.
        hexeract = (2**dimensions) * comb(dimensions+6,6)
        hexeract = "{:,}".format(hexeract)
        return hexeract

    def hepteract(dimensions):
        logging.info("Class.def 'shape.hepteract' invoked.")
        hepteract = 2**(dimensions-7)*scipy.special.binom(dimensions, 7)
        hepteract = "{:,}".format(hepteract)
        return hepteract

    def octeract(dimensions):
        logging.info("Class.def 'shape.octeract' invoked.")
        dimensions = dimensions - 8
            # For some reason, this function returns the results for 8 higher dimensions than requested. Need to research why this is.
        octeract = scipy.special.binom(dimensions+8,8) * (2**dimensions)
        octeract = "{:,}".format(octeract)
        return octeract

    def enneract(dimensions):
        logging.info("Class.def 'shape.enneract' invoked.")
        dimensions = dimensions - 9
            # For some reason, this function returns the results for 9 higher dimensions than requested. Need to research why this is.
        enneract = scipy.special.binom(dimensions+9,9) * (2**dimensions)
        enneract = "{:,}".format(enneract)
        return enneract

    def dekeract(dimensions):
        logging.info("Class.def 'shape.dekeract' invoked.")
        dimensions = dimensions - 10
            # For some reason, this function returns the results for 9 higher dimensions than requested. Need to research why this is.
        dekeract = scipy.special.binom(dimensions+9,9) * (2**dimensions)
        dekeract = "{:,}".format(dekeract)
        return dekeract

    def hypercube(dimensions):
        logging.info("Class.def 'shape.hypercube' invoked.")
        originalDimensions = dimensions
        dimensions = dimensions - dimensions
            # For some reason, this function returns the results for 9 higher dimensions than requested. Need to research why this is.
        hypercube = scipy.special.binom(dimensions+originalDimensions-1,originalDimensions-1) * (2**dimensions)
        hypercube = "{:,}".format(hypercube)
        return hypercube
     
class matrices():
    def projectionMatrix(dimensions):
        logging.info("Class.def 'matrices.projectionMatrix' invoked.")
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
        logging.info("Class.def 'matrices.verticesMatrix' invoked.")
        if dimensions == 0:
            verticesMatrix = np.array([0, 0, 0])
                # Ensures that the return value is a 3D coordinate for plotting. 
            rows = 1
            return verticesMatrix, rows
        else:
            verticesMatrix = magnitude * (2*((np.arange(2**dimensions)[:,None] & (1 << np.arange(dimensions))) > 0) - 1)
            logging.info(("Orginal VerticesMatrix Generated:", verticesMatrix))
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
                logging.error("Number of vertices doesn't match the design intent of the requested shape.")
            else:
                logging.info("Vertices successfully generated:")
                return verticesMatrix, rows
    
    def rotationMatrix(rotationAxis, theta):
        logging.info("Class.def 'matrices.rotationMatrix' invoked.")
        
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
    def flattenHypercube(dims, mag, vertices, originalVertices):
        logging.info("Class.def 'flatten.flattenHypercube' invoked.")
        
        myProjection = matrices.projectionMatrix(dims)
        logging.info(("Projection Matrix:", myProjection))
        
        myProjectionRows, myProjectionCols = np.shape(myProjection)
        myVerticesRows, myVerticesAxes = originalVertices[0].shape
        
        if myProjectionCols != myVerticesAxes:
            logging.error("# of Columns in Projection Matrix must be the same as the # of axes in each coordinate set from the Vertices Matrix.")
            logging.info("# of Cols in Projection Matrix:", myProjectionCols)
            logging.info("# of axes per Coordinate Set:", myVerticesAxes)
        else:
            logging.info("No errors preventing flattening the verticesMatrix.")
            
            myFlattenedHypercube = np.empty((myVerticesRows, myProjectionRows))
            logging.info("Shape of Empty Flattened Matrix")
            logging.info(myFlattenedHypercube.shape)
            logging.info("Empty Flattened Matrix:")
            logging.info(myFlattenedHypercube)

            row = 0
            while row <= vertices - 1:
                myFlattenedHypercube[row] = np.matmul(myProjection, originalVertices[0][row])
                logging.info((row, myFlattenedHypercube[row]))
                row += 1
            
            logging.info("Flattened Hypercube Coordinates Generated.")
            return myFlattenedHypercube


class MyMainWindow(MainWindow):
    logging.info("Class 'MyMainWindow' has been called.")   
    def __init__(self, parent=None, show=True):
        logging.info("def '__init__' has been called.")
        QtWidgets.QMainWindow.__init__(self, parent)
      
        ### CREATE FRAME ###
        logging.info("Creating frame.")
        self.frame = QtWidgets.QFrame()
        vlayout = QtWidgets.QVBoxLayout()
        
        ### PYVISTA INTERACTOR OBJECT ###
        logging.info("Creating PyVista Interactor Object.")
        self.plotter = QtInteractor(self.frame)
        vlayout.addWidget(self.plotter.interactor)
        self.signal_close.connect(self.plotter.close)

        self.frame.setLayout(vlayout)
        self.setCentralWidget(self.frame)
        
        ### MENU ###
        logging.info("Creating main menu.")
        mainMenu = self.menuBar()
        logging.info("Creating 'File' submenu.")
        fileMenu = mainMenu.addMenu('File')
        exitButton = QtWidgets.QAction('Exit', self)
        exitButton.setShortcut('Ctrl+X')
        exitButton.triggered.connect(self.close)
        fileMenu.addAction(exitButton)
        
        ### MESH SUBMENU ###
        logging.info("Creating 'Mesh' submenu.")
        meshMenu = mainMenu.addMenu('Mesh')
        
        ## Sphere ##
        logging.info("Creating 'Sphere' item in Mesh submenu.")
        self.add_sphere_action = QtWidgets.QAction('Add Sphere', self)
        self.add_sphere_action.triggered.connect(self.clear_plotter)
        self.add_sphere_action.triggered.connect(self.add_sphere)
        meshMenu.addAction(self.add_sphere_action)
        
        ## Point ##
        logging.info("Creating 'Point' item in Mesh submenu.")
        self.add_point_action = QtWidgets.QAction('Add Point', self)
        self.add_point_action.triggered.connect(self.clear_plotter)
        self.add_point_action.triggered.connect(self.add_point)
        meshMenu.addAction(self.add_point_action)
        
        ## Line ##
        logging.info("Creating 'Line' item in Mesh submenu.")
        self.add_line_action = QtWidgets.QAction('Add Line', self)
        self.add_line_action.triggered.connect(self.clear_plotter)
        self.add_line_action.triggered.connect(self.add_line)
        meshMenu.addAction(self.add_line_action)

        ## Plane ##
        logging.info("Creating 'Plane' item in Mesh submenu.")
        self.add_plane_action = QtWidgets.QAction('Add Plane', self)
        self.add_plane_action.triggered.connect(self.clear_plotter)
        self.add_plane_action.triggered.connect(self.add_plane)
        meshMenu.addAction(self.add_plane_action)

        ## Cube ##
        logging.info("Creating 'Cube' item in Mesh submenu.")
        self.add_cube_action = QtWidgets.QAction('Add Cube', self)
        self.add_cube_action.triggered.connect(self.clear_plotter)
        self.add_cube_action.triggered.connect(self.add_cube)
        meshMenu.addAction(self.add_cube_action)
        
        ## Hypercube of n=4 Dimensions ##
        logging.info("Creating 'Hypercube' 4D item in Mesh submenu.")
        self.add_hypercube_action = QtWidgets.QAction('Add Hypercube', self)
        self.add_hypercube_action.triggered.connect(self.clear_plotter)
        self.add_hypercube_action.triggered.connect(self.add_hypercube)
        meshMenu.addAction(self.add_hypercube_action)
        
        if show:
            self.show()

    def clear_plotter(self):
        logging.info("Clearing the plotter and resetting isometric view.")
        self.plotter.view_isometric()
        self.plotter.clear()

    def add_sphere(self):
        logging.info("User has requested a sphere in the plotter.")
        sphere = pv.Sphere()
        legend = "3D: A Sphere"
        self.plotter.add_mesh(sphere, show_edges=True, label=legend)
        self.plotter.add_axes(line_width=3, labels_off=False)
        self.plotter.reset_camera()
      
    def add_point(self):
        logging.info("User has requested a point in the plotter.")
        myNewVertices = matrices.verticesMatrix(0, 1)
        numberOfPoints = shape.points(0)
        myLabels = ()
        for rows in range(numberOfPoints):
            myLabels = np.append(myLabels, rows + 1)
        legend = "OD: A Point"
        self.plotter.add_mesh(pv.PolyData(myNewVertices[0]), smooth_shading=True, label=legend, color='white', show_edges=True)
        self.plotter.add_point_labels(myNewVertices[0], myLabels, font_size=10, point_color='red', point_size=10, render_points_as_spheres=True, always_visible=True, fill_shape=False)
        self.plotter.add_axes(line_width=3, labels_off=False)
        self.plotter.add_legend(face=None)

    def add_line(self):
        logging.info("User has requested a line in the plotter.")
        myNewVertices = matrices.verticesMatrix(1, 1)
        numberOfPoints = shape.points(1)
        myLabels = ()
        for rows in range(numberOfPoints):
            myLabels = np.append(myLabels, rows + 1)
        legend = "1D: A Line"
        self.plotter.add_mesh(pv.PolyData(myNewVertices[0]), smooth_shading=True, label=legend, color='white', show_edges=True)
        self.plotter.add_point_labels(myNewVertices[0], myLabels, font_size=10, point_color='red', point_size=10, render_points_as_spheres=True, always_visible=True, fill_shape=False)
        self.plotter.add_axes(line_width=3, labels_off=False)
        self.plotter.add_legend(face=None)
        i = 0
        lines = np.array((myNewVertices[0][i], myNewVertices[0][i+1]))
        self.plotter.add_lines(lines)

    def add_plane(self):
        logging.info("User has requested a plane in the plotter.")
        myNewVertices = matrices.verticesMatrix(2, 1)
        numberOfPoints = shape.points(2)
        myLabels = ()
        for rows in range(numberOfPoints):
            myLabels = np.append(myLabels, rows + 1)
        legend = "2D: A Plane"
        self.plotter.add_mesh(pv.PolyData(myNewVertices[0]), smooth_shading=True, label=legend, color='white', show_edges=True)
        self.plotter.add_point_labels(myNewVertices[0], myLabels, font_size=10, point_color='red', point_size=10, render_points_as_spheres=True, always_visible=True, fill_shape=False)
        self.plotter.add_axes(line_width=3, labels_off=False)
        self.plotter.add_legend(face=None)
        i = 0
        while i < 4: # (dim * 2)
            lines = np.array((myNewVertices[0][i], myNewVertices[0][(i + 2) % 2])) # ((i + dim) % dim)
            self.plotter.add_lines(lines)
            i += 1 # (dim / 2)
        i = 0
        while i < 4: # (dim * 2)
            lines2 = np.array((myNewVertices[0][i], myNewVertices[0][(i + 1) % 4])) # (i + (dim / 2)) % (dim / 2))
            self.plotter.add_lines(lines2)
            i += 2 # (dim)
        # Draw four lines connecting the plane.
      
    def add_cube(self):
        logging.info("User has requested a cube in the plotter.")
        myNewVertices = matrices.verticesMatrix(3, 1)
        numberOfPoints = shape.points(3)
        myLabels = ()
        for rows in range(numberOfPoints):
            myLabels = np.append(myLabels, rows + 1)
        legend = "3D: A Cube"
        self.plotter.add_mesh(pv.PolyData(myNewVertices[0]), smooth_shading=True, label=legend, color='white', show_edges=True)
        self.plotter.add_point_labels(myNewVertices[0], myLabels, font_size=10, point_color='red', point_size=10, render_points_as_spheres=True, always_visible=True, fill_shape=False)
        self.plotter.add_axes(line_width=3, labels_off=False)
        self.plotter.add_legend(face=None)
        i = 0
        while i < 4: # (dim - 1)
            lines = np.array((myNewVertices[0][i], myNewVertices[0][(i+2)%4])) # (i + (dim - 1) % (dim + 1))
            self.plotter.add_lines(lines)
            lines2 = np.array((myNewVertices[0][i+4], myNewVertices[0][((i+2)%4)+4])) # ((i + (dim - 1)) % (dim + 1) + (dim + 1)) 
            self.plotter.add_lines(lines2)
            lines3 = np.array((myNewVertices[0][i], myNewVertices[0][(i+4)])) # (i + (dim + 1))
            self.plotter.add_lines(lines3)
            i += 1 # (dim / 3)
        i = 0
        while i < 4: # (dim - 1)
            lines4 = np.array((myNewVertices[0][i], myNewVertices[0][(i+1)%4])) # (i + (dim - 3) % (dim + 1))
            self.plotter.add_lines(lines4)
            lines5 = np.array((myNewVertices[0][i+4], myNewVertices[0][((i+1)%4)+4])) # ((i + (dims - 3) % (dim + 1)) + (dim + 1))
            self.plotter.add_lines(lines5)
            i += 2 # (dim - 1)
            # Draw lines connecting two planes to each other.

    def add_hypercube(self):
        logging.info("User has requested a hypercube - 4D in the plotter.")
        myNewVertices = matrices.verticesMatrix(4, 1)
        numberOfPoints = shape.points(4)

        myFlattenedVertices = flatten.flattenHypercube(4, 1, numberOfPoints, myNewVertices)
        logging.info(("Returned flattened vertices in def 'add_hypercube':", myFlattenedVertices))

        myLabels = ()
        for rows in range(numberOfPoints):
            myLabels = np.append(myLabels, rows + 1)
        legend = "4D: A Hypercube of n=4 Dimensions"
        logging.info(("Hypercube Point Labels:", myLabels))
        self.plotter.add_mesh(pv.PolyData(myFlattenedVertices), smooth_shading=True, label=legend, color='white', show_edges=True)
        #self.plotter.add_point_labels(myFlattenedVertices[0], myLabels, font_size=8, point_color='red', point_size=8, render_points_as_spheres=True, always_visible=True, fill_shape=False)
        # The second half of the points generated automatically are identicle to the first set of points generated, so instead of getting 16 unique points, I am
        # getting 8 points, twice.
        self.plotter.add_axes(line_width=3, labels_off=False)
        self.plotter.add_legend(face=None)
       
### RUN ###
if __name__ == "__main__":
    logging.info("executing 'if __name__ == __main__' to start app.")
    app = QtWidgets.QApplication(sys.argv)
    window = MyMainWindow()
    sys.exit(app.exec_())
   
else:
    logging.CRITICAL("This module is not __main__ and should not be directly executed.")
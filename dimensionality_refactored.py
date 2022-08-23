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

os.environ["QT_API"] = "pyqt5"

def loggingAndWarnings():
   ### LOGGING ###
   lvl = logging.DEBUG
   fmt = '[%(levelname)s] %(asctime)s :: %(message)s'
   logging.basicConfig(level=lvl, format=fmt)

   ### WARNINGS HANDLERS ###
   warnings.filterwarnings('ignore', '.*force_float.*')
   warnings.filterwarnings('ignore', '.*auto_close*')

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

class MyMainWindow(MainWindow):   
   def __init__(self, parent=None, show=True):
      QtWidgets.QMainWindow.__init__(self, parent)
      
      ### CREATE FRAME ###
      self.frame = QtWidgets.QFrame()
      vlayout = QtWidgets.QVBoxLayout()
      
      ### PYVISTA INTERACTOR OBJECT ###
      self.plotter = QtInteractor(self.frame)
      vlayout.addWidget(self.plotter.interactor)
      self.signal_close.connect(self.plotter.close)

      self.frame.setLayout(vlayout)
      self.setCentralWidget(self.frame)
      
      ### MENU ###
      mainMenu = self.menuBar()
      fileMenu = mainMenu.addMenu('File')
      exitButton = QtWidgets.QAction('Exit', self)
      exitButton.setShortcut('Ctrl+X')
      exitButton.triggered.connect(self.close)
      fileMenu.addAction(exitButton)
      
      ### MESH SUBMENU ###
      meshMenu = mainMenu.addMenu('Mesh')
      
      ## Sphere ##
      self.add_sphere_action = QtWidgets.QAction('Add Sphere', self)
      self.add_sphere_action.triggered.connect(self.clear_plotter)
      self.add_sphere_action.triggered.connect(self.add_sphere)
      meshMenu.addAction(self.add_sphere_action)
      
      ## Point ##
      self.add_point_action = QtWidgets.QAction('Add Point', self)
      self.add_point_action.triggered.connect(self.clear_plotter)
      self.add_point_action.triggered.connect(self.add_point)
      meshMenu.addAction(self.add_point_action)
      
      ## Line ##
      self.add_line_action = QtWidgets.QAction('Add Line', self)
      self.add_line_action.triggered.connect(self.clear_plotter)
      self.add_line_action.triggered.connect(self.add_line)
      meshMenu.addAction(self.add_line_action)

      ## Plane ##
      self.add_plane_action = QtWidgets.QAction('Add Plane', self)
      self.add_plane_action.triggered.connect(self.clear_plotter)
      self.add_plane_action.triggered.connect(self.add_plane)
      meshMenu.addAction(self.add_plane_action)

      ## Cube ##
      self.add_cube_action = QtWidgets.QAction('Add Cube', self)
      self.add_cube_action.triggered.connect(self.clear_plotter)
      self.add_cube_action.triggered.connect(self.add_cube)
      meshMenu.addAction(self.add_cube_action)
      
      if show:
         self.show()

   def clear_plotter(self):
      self.plotter.view_isometric()
      self.plotter.clear()

   def add_sphere(self):
      sphere = pv.Sphere()
      legend = "3D: A Sphere"
      self.plotter.add_mesh(sphere, show_edges=True, label=legend)
      self.plotter.add_axes(line_width=3, labels_off=False)
      self.plotter.reset_camera()
      
   def add_point(self):
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
      
### RUN ###
if __name__ == "__main__":
   logWarn = threading.Thread(target=loggingAndWarnings(), daemon=True)
   logWarn.start()
   
   app = QtWidgets.QApplication(sys.argv)
   window = MyMainWindow()
   sys.exit(app.exec_())
   
else:
   logging.debug("main() failed to initialize properly.")
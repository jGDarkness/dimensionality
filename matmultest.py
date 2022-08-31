import numpy as np

dims = 4
magnitude = 1

# Projection Matrix Generator
projectionMatrix = np.eye(dims, dims, dtype=np.int32)
print ('\n', projectionMatrix, '\n', projectionMatrix.shape)

# Vertices Generator
if dims == 0:
   verticesMatrix = np.array([0, 0, 0])
elif dims == 1:
   verticesMatrix = magnitude * (2*((np.arange(2**1)[:,None] & (1 << np.arange(1))) > 0) - 1)
   remainingAxes = np.array([0,0])
   verticesMatrix = np.insert(verticesMatrix,1,remainingAxes,axis=1)
   verticesMatrix = np.insert(verticesMatrix,2,remainingAxes,axis=1)
elif dims == 2:
   verticesMatrix = magnitude * (2*((np.arange(2**2)[:,None] & (1 << np.arange(2))) > 0) - 1)
   remainingAxes = np.array([0,0,0,0])
   verticesMatrix = np.insert(verticesMatrix,2,remainingAxes,axis=1)          
else:
   verticesMatrix = magnitude * (2*((np.arange(2**dims)[:,None] & (1 << np.arange(dims))) > 0) - 1)

print ('\n', verticesMatrix, '\n', verticesMatrix.shape)

if dims >= 4:
   rows = 0
   while rows < 16:
      flattenedVertices = np.matmul(projectionMatrix, verticesMatrix[0][rows])
      rows += 1
      
   
   print ('\n', flattenedVertices, '\n', flattenedVertices.shape)



print ('\n')
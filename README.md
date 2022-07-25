# dimensionality
 A python study of dimensionality

## Purpose
 I am fully aware that there are many and varied software, language, and mathematical models fully available on the market to enable the kind of study I am doing, however, part of my purpose is to gain a greater understanding of object oriented programming, Python, and interactive software design.

## Dimensions
 My motivation for this particular field of study, is based on reading the novels in the "Three Body" universe written by Cixin Liu. And while his novels deal with dimensionality in the realm of standard and quantum physics, including string theory, I am primarily interested in the ability to "unfold" higher dimensions into lower dimensions, i.e., rendering a 3D model of an 10D shape. As a result, I am limiting myself to the Euclidean geometrical shapes of a point (0D), a line, a face, a cube, a hypercube or tesseract, and so on to a dekeract (10D). I hope to be able to develop models of higher dimensionality (>10D), however, reverse engineering the number sequence formulae for dimensions 10 and higher may be beyond my skillset.

## Source Materials
 I have made liberal use of the Online Encyclopedia of Integer Sequences in order to put this study together, and other principles of Euclidean Geometry and the follow on studies made by important mathematicians.

## Data Tables and Formulae
 This datatable and accompanying formulae are how I can validate the results of my calculations in various parts of the program.

The number E_sub(m,n) of m-dimensional faces of a n-dimensional hypercube
                         0face  1face  2face  3face  4face 5face 6face 7face 8face 9face 10face
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

The number of "0-face"  Points (Vertex)                in each Array[1]  is 2^n                                 
The number of "1-face"  Lines Segments (Edge)          in each Array[2]  is n*2^(n-1)                           
The number of "2-face"  Squares (Face)                 in each Array[3]  is where n = n-1, let n*(n+1)*2^(n-2)    offset 1 zeroes
The number of "3-face"  Cubes (Cell)                   in each Array[4]  is binomial[n, 3]*2^(n-3)             
The number of "4-face"  Tesseracts (hypercube)         in each Array[5]  is 2^(n-4)*C(n,4)                      
The number of "5-face"  Penteracts (5D hypercubes)     in each Array[6]  is 2^(n-5)*binomial(n,5)               
The number of "6-face"  Hexeracts (n+6D hypercubes)    in each Array[7]  is where n = n-5, let 2^n*C(n+6,6)       offset 6 zeroes
The number of "7-face"  Hepteracts (n+7D hypercubes)   in each Array[8]  is 2^(n-7)*binomial(n,7)               
The number of "8-face"  Octeracts (n+8D hypercubes)    in each Array[9]  is binomial(n+8,8) * 2^n                 offset 8 zeroes
The number of "9-face"  Enneracts (n+9D hypercubes)    in each Array[10] is binomial(n+9,9) * 2^n                 offset 9 zeroes
The number of "10-face" Dekeract (n+10D hypercubes)    in each Array[11] is binomial(0+n-1,n-1) * 2^0             offset 10 or more zeroes based on dimensionality
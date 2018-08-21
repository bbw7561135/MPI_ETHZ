FILES
cpp_ex3.cpp cpp_ex3.h

FUNCTION
The function of the code is not important.  However, so you have a basic
understanding of what is going on, it is merely taking a House object made of
Triangles and Rectangles and moving it back and forth along a line in 3D
space creating a "movie".  There are no graphics and no actual movie is
created but merely an STL vector of House objects.

EXPECTED OUTPUT
    Front: [(0,0,0),(0,1,0),(1,1,0),(1,0,0)],[(0,1,0),(1,1,0),(0.5,1.5,0)]
    Back:  [(0,0,1),(0,1,1),(1,1,1),(1,0,1)],[(0,1,1),(1,1,1),(0.5,1.5,1)]
    Left:
    [(0,0,0),(0,1,0),(0,1,1),(0,0,1)],[(0,1,0),(0.5,1.5,0),(0.5,1.5,1),(0,0,1)]
    Right:
    [(1,1,0),(1,0,0),(1,0,1),(1,1,1)],[(1,1,0),(0.5,1.5,0),(0.5,1.5,1),(1,1,1)]


    The following should be the same
    Front: [(0,0,0),(0,1,0),(1,1,0),(1,0,0)],[(0,1,0),(1,1,0),(0.5,1.5,0)]
    Back:  [(0,0,1),(0,1,1),(1,1,1),(1,0,1)],[(0,1,1),(1,1,1),(0.5,1.5,1)]
    Left:
    [(0,0,0),(0,1,0),(0,1,1),(0,0,1)],[(0,1,0),(0.5,1.5,0),(0.5,1.5,1),(0,0,1)]
    Right:
    [(1,1,0),(1,0,0),(1,0,1),(1,1,1)],[(1,1,0),(0.5,1.5,0),(0.5,1.5,1),(1,1,1)]

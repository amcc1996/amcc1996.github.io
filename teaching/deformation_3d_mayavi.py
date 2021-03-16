# Master in Mechanical Engineering
# Solid Mechanics
#                                   Faculty of Engineering of University of Porto, 2019/2020
# ==========================================================================================
# Deformed configuration in Problem 1 on Large Deformations using Mayavi.
#
#                                           A. M. Couto Carneiro <amcc@fe.up.pt>, @FEUP/CM2S
# ==========================================================================================
# Import modules
# --------------
# Array computation.
import numpy as np
# Mayavi for nice 3D plotting
from mayavi import mlab

# Functions
# ---------
def deformation(X1, X2, X3):
    x1 = X1 - A * X3
    x2 = X2 - A * X3
    x3 = -A * X1 + A * X2 + X3
    return x1, x2, x3
# ==========================================================================================
A = 0.5
N = 2
a = 0.2

side = np.linspace(-1, 1, num=N, endpoint=True)

colorref = (228.0/256.0, 26.0/256.0, 28.0/256.0)
opacref = 0.6
colordef = (55.0/256.0, 126.0/256.0, 184.0/256.0)
opacdef = 1.0
ref = "surface"

fig = mlab.figure(size=(1000, 1000))

# Bottom face
X1, X2 = np.meshgrid(side, side, indexing='xy')
X3 = np.full((N, N), -1)
mlab.mesh(X1, X2, X3, color=colorref, opacity=opacref, representation=ref)
x1, x2, x3 = np.vectorize(deformation)(X1, X2, X3)
mlab.mesh(x1, x2, x3, color=colordef, opacity=opacdef, representation=ref)

# Top face
X1, X2 = np.meshgrid(side, side, indexing='xy')
X3 = np.full((N, N), 1)
mlab.mesh(X1, X2, X3, color=colorref, opacity=opacref, representation=ref)
x1, x2, x3 = np.vectorize(deformation)(X1, X2, X3)
mlab.mesh(x1, x2, x3, color=colordef, opacity=opacdef, representation=ref)

# Back face
X2, X3 = np.meshgrid(side, side, indexing='xy')
X1 = np.full((N, N), -1)
mlab.mesh(X1, X2, X3, color=colorref, opacity=opacref, representation=ref)
x1, x2, x3 = np.vectorize(deformation)(X1, X2, X3)
mlab.mesh(x1, x2, x3, color=colordef, opacity=opacdef, representation=ref)

# Front face
X2, X3 = np.meshgrid(side, side, indexing='xy')
X1 = np.full((N, N), 1)
mlab.mesh(X1, X2, X3, color=colorref, opacity=opacref, representation=ref)
x1, x2, x3 = np.vectorize(deformation)(X1, X2, X3)
mlab.mesh(x1, x2, x3, color=colordef, opacity=opacdef, representation=ref)

# Left face
X3, X1 = np.meshgrid(side, side, indexing='xy')
X2 = np.full((N, N), -1)
mlab.mesh(X1, X2, X3, color=colorref, opacity=opacref, representation=ref)
x1, x2, x3 = np.vectorize(deformation)(X1, X2, X3)
mlab.mesh(x1, x2, x3, color=colordef, opacity=opacdef, representation=ref)

# Right face
X3, X1 = np.meshgrid(side, side, indexing='xy')
X2 = np.full((N, N), 1)
mlab.mesh(X1, X2, X3, color=colorref, opacity=opacref, representation=ref)
x1, x2, x3 = np.vectorize(deformation)(X1, X2, X3)
mlab.mesh(x1, x2, x3, color=colordef, opacity=opacdef, representation=ref)

mlab.show()

# Master in Mechanical Engineering
# Solid Mechanics
#                                   Faculty of Engineering of University of Porto, 2019/2020
# ==========================================================================================
# Graphical illutration of shear strains using the displacement field associated with a
# unixial tensile loading. There is no shear strain associated with the direction of the
# loading and the in-plane orthogonal direction, but if we analyse the deformation of
# segemnts at +-45º the shear strain is evidenced.
#
#                                           A. M. Couto Carneiro <amcc@fe.up.pt>, @FEUP/CM2S
# ==========================================================================================
# Import modules
# --------------
# Matplotlib for plotting and producing the animation
from mpl_toolkits.mplot3d import Axes3D, proj3d
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch
from matplotlib.patches import Circle, Wedge, Polygon
from matplotlib.collections import PatchCollection
from matplotlib import animation
from matplotlib import rc
# Simple array computational toolbox
import numpy as np
from numpy.linalg import norm
# Basic mathematical toolbox
from math import cos, sin, pi
# Time handling module
import time
# ==========================================================================================
def disp(x, y, step):
    """Displacement function.
    """
    poisson = 0.3
    exx = 1 * par[step]
    eyy =  - poisson * exx
    u = exx * x
    v = eyy * y
    return u, v

def rotate(x, y):
    """Rotate a (x, y) point relative to the origin of the coordinate system.
    """
    xr = x * np.cos(angle) - y * np.sin(angle)
    yr = x * np.sin(angle) + y * np.cos(angle)
    return xr, yr

# Plot settings
# -------------
rc("font", family="serif")
rc("font", size=22)
rc('text', usetex=False)
plt.style.use('dark_background')
# ==========================================================================================
# Grid definition.
angle = np.pi/4.0 * 1
nx = 20
ny = 20
x = np.linspace(-1, 1, num=nx, endpoint=True)
y = np.linspace(-1, 1, num=ny, endpoint=True)
xg, yg = np.meshgrid(x, y, indexing='xy')
xg, yg = np.vectorize(rotate)(xg, yg)

# deformation parameterization
par = np.linspace(0, 1, num=200, endpoint=True)

# Create figure and axis.
fig = plt.figure(num='Shear strain demonstration', figsize=(8, 8), constrained_layout=True)
ax = fig.add_subplot(111)

# Save animation.
save = False

def init():
    return ax

def animate(i):
    plt.cla()
    u, v = np.vectorize(disp)(xg, yg, i)
    ax.pcolormesh(xg + u, yg + v, np.zeros((ny, nx)), edgecolor='white', cmap="gist_gray")
    ax.scatter(xg + u, yg + v, s=80, c='cornflowerblue', marker="o")
    ax.axis("off")
    ax.set_aspect('equal', adjustable='box')
    ax.set_xlim(-0.25, 0.25)
    ax.set_ylim(-0.25, 0.25)

# Create the animation.
anim = animation.FuncAnimation(fig, animate, init_func=init,
                               frames=len(par), interval=50, blit=False, repeat=False)

if save:
    Writer = animation.writers['ffmpeg']
    writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=1800)
    anim.save('shearstrain45.mp4', writer=writer)

plt.show()

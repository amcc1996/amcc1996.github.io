# Master in Mechanical Engineering
# Solid Mechanics
#                                   Faculty of Engineering of University of Porto, 2019/2020
# ==========================================================================================
# Principal stresses and principal directions. This script generates and animation
# of a rotating normal vector and the respective stress tensor. As the normal approaches a
# principal direction, and ultimately matches the principal direction the stress vector
# becomes colinear with the normal vector.
#
#                                           A. M. Couto Carneiro <amcc@fe.up.pt>, @FEUP/CM2S
# ==========================================================================================
# Import modules
# --------------
# Matplotlib for plotting and producing the animation
from mpl_toolkits.mplot3d import Axes3D, proj3d
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch
from matplotlib import animation
from matplotlib import rc
# Simple array computational toolbox
import numpy as np
from numpy.linalg import norm
# Basic mathematical toolbox
from math import cos, sin, pi
# ==========================================================================================
# Plot settings
# -------------
rc("font",family="serif")
rc("font",size=18)
rc('text', usetex=True)
plt.style.use('dark_background')
# ================================================================================== Arrow3D
class Arrow3D(FancyArrowPatch):
    """3D Arrow object.
    """
    def __init__(self, xs, ys, zs, *args, **kwargs):
        FancyArrowPatch.__init__(self, (0,0), (0,0), *args, **kwargs)
        self._verts3d = xs, ys, zs

    def draw(self, renderer):
        xs3d, ys3d, zs3d = self._verts3d
        xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, renderer.M)
        self.set_positions((xs[0],ys[0]),(xs[1],ys[1]))
        FancyArrowPatch.draw(self, renderer)
# =============================================================================== drawVector
def drawVector(vec, color, lines=True, head = '-|>', linewidth=2):
    """Simplified function to draw the vector.
    """
    arrow = Arrow3D([0, vec[0]], [0, vec[1]],[0, vec[2]], mutation_scale=20,
                    lw=linewidth, arrowstyle=head, color=color, zorder=10)

    if lines:
        ax.plot([vec[0], vec[0]], [vec[1], vec[1]], zs=[vec[2], 0], color=color,
                 linewidth=1)
        ax.plot([vec[0], 0], [vec[1], vec[1]], zs=[0, 0], color=color, linewidth=1)
        ax.plot([vec[0], vec[0]], [vec[1], 0], zs=[0, 0], color=color, linewidth=1)
    return arrow
# ==========================================================================================
# Principal direction to plot.
pDir = 3

# Cauchy stress tensor.
sigma = np.array([[0, 1, 0],
                  [1, 2, 2],
                  [0, 2, 0]])

# Known principal directions.
n1f = np.array([0.243, 0.839, 0.486])
n2f = np.array([0.894, 0.000, -0.447])
n3f = np.array([0.375, -0.544, 0.750])

# Parameter spaces.
par = np.linspace(0, 1, num = 200, endpoint= True)

# Create figure and axis.
fig = plt.figure(num='Tensões e direções principais', figsize = (9,7))
ax = fig.gca(projection='3d')
plt.tight_layout()

# Scale constant to set the length of the coordinate system
scale = 3

# Save animations
save = False
# ------------------------------------------------------------------------------------------
def init():
    return ax
# ------------------------------------------------------------------------------------------
def animate(i):
    # Clear axis.
    ax.cla()

    # Draw origin of coordinate system.
    ax.scatter([0],[0],[0],color="white",s=50)

    # Setup coordinate system.
    x = Arrow3D([0, scale], [0, 0], [0, 0], mutation_scale=20, lw=1, arrowstyle="-|>",
                color='white', zorder=1)
    ax.plot([0, -scale], [0 ,0], zs=[0, 0], color='white', linewidth=1, linestyle='--')

    y = Arrow3D([0, 0], [0, scale], [0, 0], mutation_scale=20, lw=1, arrowstyle="-|>",
                color='white', zorder=1)
    ax.plot([0, 0], [0 ,-scale], zs=[0, 0], color='white', linewidth=1, linestyle='--')

    z = Arrow3D([0, 0], [0, 0], [0, scale], mutation_scale=20, lw=1, arrowstyle="-|>",
                color='white', zorder=1)
    ax.plot([0, 0], [0 ,0], zs=[0, -scale], color='white', linewidth=1, linestyle='--')

    # Axis labels.
    ax.text(scale, 0, -0.2, r'$x$', None)
    ax.text(0, scale, -0.2, r'$y$', None)
    ax.text(0, 0, scale, r'$z$', None)

    # Initial view.
    ax.view_init(elev=45, azim=40)

    # Draw axis.
    ax.add_artist(x)
    ax.add_artist(y)
    ax.add_artist(z)

    # Do not dar the default 3D box.
    ax.set_axis_off()

    # Compute the normal vector of the current frame.
    # The components of each vector are parameterized such that when t=1, this is, in the
    # last frame, the normal vector equals the respective principal direction.
    n1 = np.array([n1f[0] * cos(par[i] * 2 * pi),
                  n1f[1] + sin(par[i] * pi) / 4,
                  n1f[2] + (par[i] - 1) ** 2])

    n2 = np.array([n2f[0] * cos(par[i] * 2 * pi),
                  n2f[1] + sin(par[i] * pi) / 4,
                  n2f[2] + (par[i] - 1) ** 2])

    n3 = np.array([n3f[0] * cos(par[i] * 2 * pi),
                  n3f[1] + sin(par[i] * pi) / 4,
                  n3f[2] + (par[i] - 1) ** 2])


    # Normalize the normal vector.
    n1 = n1 / norm(n1)
    n2 = n2 / norm(n2)
    n3 = n3 / norm(n3)

    # Compute the current stress vector.
    T1 = np.matmul(sigma, n1)
    T2 = np.matmul(sigma, n2)
    T3 = np.matmul(sigma, n3)

    # Set up the arrows for the current normal and stress vector.
    if pDir == 1:
        narrow1 = drawVector(n1, 'deepskyblue')
        Tarrow1 = drawVector(T1, 'red')

        ax.add_artist(narrow1)
        ax.add_artist(Tarrow1)

        ax.text(n1[0], n1[1], n1[2], '$\\vec{n}$', None,
                bbox=dict(facecolor='white', alpha=0.2, edgecolor='white'))
        ax.text(T1[0], T1[1], T1[2], '$\\vec{T}$', None,
                bbox=dict(facecolor='white', alpha=0.2, edgecolor='white'))
    elif pDir == 2:
        narrow2 = drawVector(n2, 'deepskyblue')
        Tarrow2 = drawVector(T2, 'red')

        ax.add_artist(narrow2)
        ax.add_artist(Tarrow2)

        ax.text(n2[0], n2[1], n2[2], '$\\vec{n}$', None,
                bbox=dict(facecolor='white', alpha=0.2, edgecolor='white'))
        ax.text(T2[0], T2[1], T2[2], '$\\vec{T}$', None,
                bbox=dict(facecolor='white', alpha=0.2, edgecolor='white'))
    elif pDir ==3:
        narrow3 = drawVector(n3, 'deepskyblue')
        Tarrow3 = drawVector(T3, 'red')

        ax.add_artist(narrow3)
        ax.add_artist(Tarrow3)

        ax.text(n3[0], n3[1], n3[2], '$\\vec{n}$', None,
                bbox=dict(facecolor='white', alpha=0.2, edgecolor='white'))
        ax.text(T3[0], T3[1], T3[2], '$\\vec{T}$', None,
                bbox=dict(facecolor='white', alpha=0.2, edgecolor='white'))

    # Set axis limits.
    ax.set_xlim(-1.5, 1.5)
    ax.set_ylim(-1.5, 1.5)
    ax.set_zlim(-1.5, 1.5)

    return ax
# ------------------------------------------------------------------------------------------
# Create the animation.
anim = animation.FuncAnimation(fig, animate, init_func=init,
                               frames=len(par), interval=50, blit=False, repeat=False)

if save:
    Writer = animation.writers['ffmpeg']
    writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=1800)
    anim.save('principalStressesN3.mp4', writer=writer)


plt.show()

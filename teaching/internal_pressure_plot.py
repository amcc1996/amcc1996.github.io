# Master of Management and Industrial Engineering
# Strength of Materials
#
# Stresses in a cylinder subjected to internal pressure: thin-walled solution and Lamé
# equations.
#                                   Faculty of Engineering of University of Porto, 2019/2020
# ==========================================================================================
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
# Plot settings
# -------------
rc("font", family="serif")
rc("font", size=12)
rc('text', usetex=False)
plt.style.use('dark_background')
# ==========================================================================================
# coordinate vectors
N = 200
ratioi = 1/1.01
ratiof = 1/100
delta = np.linspace(ratioi, ratiof, num=N, endpoint=True)
re = 0.15
ri = (1 - delta) * re
alpha = np.linspace(0, 1, num=200, endpoint=True)

plotRadial = True

# Create figure and axis.
fig = plt.figure(num='Cylinder subjected to internal pressure', figsize=(8, 8), constrained_layout=True)
ax = fig.add_subplot(111)

props = dict(boxstyle='round', facecolor=None, alpha=0.0)

save = False

def init():
    return ax

def animate(i):
    plt.cla()
    r = alpha * (re - ri[i]) + ri[i]
    outer = Circle((0.25, 1.25), re, color='silver', alpha=0.6)
    inner = Circle((0.25, 1.25), ri[i], color='black')
    ax.add_artist(outer)
    ax.add_artist(inner)
    if plotRadial:
        sr = - ri[i] ** 2 / (re ** 2 - ri[i] ** 2) * (1 - re ** 2 / r ** 2)
        ax.plot(alpha, sr)
        ax.set_ylim(0, 1.5)
        ax.plot([0, 1], [1, 0], linestyle="--", color="plum", linewidth=1.5)
        ax.fill_between(alpha, 0, sr, color="deepskyblue", alpha=0.25)
        ax.set_ylim(0, 1.5)
        ax.text(0.45, 0.39, "Thin-wall solution", transform=ax.transAxes, fontsize=14,
        verticalalignment='top', bbox=props, rotation=-45)
        ax.set_xlabel(r"$\dfrac{{r - r_i}}{{r_e-r_i}}$")
        ax.set_ylabel(r"$\dfrac{{-\sigma_r}}{{p}}$")
    else:
        st = ri[i] ** 2 / (re ** 2 - ri[i] ** 2) * (1 + re ** 2 / r ** 2) * (re - ri[i]) / (2 * re)
        ax.plot(alpha, st, color="deepskyblue", linewidth=3.0)
        ax.plot([0, 1], [0.5, 0.5], linestyle="--", color="plum", linewidth=1.5)
        ax.fill_between(alpha, 0, st, color="deepskyblue", alpha=0.25)
        ax.set_ylim(0, 1.5)
        ax.text(0.62, 0.36, "Thin-wall solution", transform=ax.transAxes, fontsize=14,
        verticalalignment='top', bbox=props)
        ax.set_xlabel(r"$\dfrac{{r - r_i}}{{r_e-r_i}}$")
        ax.set_ylabel(r"$\dfrac{{\sigma_t}}{{p}} \cdot \dfrac{{r_e - r_i}}{{2r_e}}$")

    ax.set_xlim(0, 1)

    ax.set_aspect('equal', 'box')

    # place a text box in upper left in axes coords
    ax.text(0.5, 0.85, r"$\dfrac{{r_e -  r_i}}{{2r_e}} = \dfrac{{t}}{{r_e}}  = {0:=9.4}$".format(delta[i]), transform=ax.transAxes, fontsize=14,
    verticalalignment='top', bbox=props)

# Create the animation.
anim = animation.FuncAnimation(fig, animate, init_func=init,
                               frames=N, interval=50, blit=False, repeat=False)

if save:
    Writer = animation.writers['ffmpeg']
    writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=1800)
    if plotRadial:
        anim.save('internalPressureRadial.mp4', writer=writer)
    else:
        anim.save('internalPressureCirc.mp4', writer=writer)

plt.show()

# Master in Mechanical Engineering
# Solid Mechanics
#                                   Faculty of Engineering of University of Porto, 2019/2020
# ==========================================================================================
# When solving a plane strain problem using Mohr's stress circle, one shall plot
# the points (e_xx, -gamma_xy/2) and (e_yy, gamma_xy/2), which might, at first
# seem rather odd. By inspecting the equations for the normal and shear stress
# on a rotated direction closely, it can be observed that the parameterization of
# Mohr's strain circle in terms of (e, gamma) implies that the points moves
# in the clockwise direction over the circle when the axis rotates in the
# counter clockwise direction. In order to match the direction of rotation of the
# two objects, one shall plot the point (e_xx, -gamma_xy/2), in which case both the
# point and the axis rotate counter clockwise.
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
# Plot settings
# -------------
rc("font", family="serif")
rc("font", size=22)
rc('text', usetex=True)
plt.style.use('dark_background')
# ==========================================================================================
# Cauchy stress tensor for a plane stress problem.
sxx = 100
txy = 50
syy = 50

# Parameter space.
par = np.linspace(0, 31.72 * pi / 180, num=200, endpoint=True)
par2 = np.linspace(0, pi, num=1000, endpoint=True)

# Create figure and axis.
fig = plt.figure(num='Círculo de Mohr', figsize=(12, 8))
ax = fig.add_subplot(111)
plt.tight_layout()

# Set rotation direction.
ccw = 1

# Save animation.
save = False

def init():
    return ax

def animate(i):

    if i == 4:
        time.sleep(4)
    # Clear axis.
    ax.cla()

    # Compute normal and shear stress.
    sxx2 = (sxx + syy) / 2 + (sxx - syy) / 2 * \
        np.cos(2 * par[i]) + txy * np.sin(2 * par[i])
    syy2 = (sxx + syy) / 2 + (sxx - syy) / 2 * np.cos(2 *
                                                      (par[i]+pi/2)) + txy * np.sin(2 * (par[i] + pi/2))
    txy2 = - (sxx - syy) / 2 * np.sin(2 * par[i]) + txy * np.cos(2 * par[i])

    ax.set_xlim(-30, 170)
    ax.set_ylim(-80, 80)

    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()

    # Removing the default axis on all sides.
    for side in ['bottom', 'right', 'top', 'left']:
        ax.spines[side].set_visible(False)

    # Removing the axis ticks:
    plt.xticks([])  # labels;
    plt.yticks([])
    ax.xaxis.set_ticks_position('none')  # tick markers.
    ax.yaxis.set_ticks_position('none')

    # Get width and height of axes object to compute matching arrowhead length and width.
    dps = fig.dpi_scale_trans.inverted()
    bbox = ax.get_window_extent().transformed(dps)
    width, height = bbox.width, bbox.height

    # Manual arrowhead width and length.
    hw = 1./40.*(ymax-ymin)
    hl = 1./40.*(xmax-xmin)
    lw = 1.  # axis line width
    ohg = 0.3  # arrow overhang

    # Compute matching arrowhead length and width.
    yhw = hw/(ymax-ymin)*(xmax-xmin) * height/width
    yhl = hl/(xmax-xmin)*(ymax-ymin) * width/height

    # Draw sigma and tau axis.
    ax.arrow(xmin, 0, xmax-xmin, 0., fc='w', ec='w', lw=lw, head_width=hw, head_length=hl,
             overhang=ohg, length_includes_head=True, clip_on=False)
    #
    ax.arrow(0, ymin, 0., ymax-ymin, fc='w', ec='w', lw=lw, head_width=hw,
             head_length=hl, overhang=ohg, length_includes_head=True, clip_on=False)

    ax.set_aspect('equal', 'box')

    # Add labels.
    ax.text(165, -10, r'$\varepsilon$', fontsize=22)
    ax.text(-10, 75, r'$\frac{1}{2}\gamma$', fontsize=22)

    if ccw == 1:
        first = r"$(\varepsilon_{x'x'},\,-\frac{1}{2}\gamma_{x'y'})$"
        second = r"$(\varepsilon_{y'y'},\,\frac{1}{2}\gamma_{x'y'})$"
    else:
        first = r"$(\varepsilon_{x'x'},\,\frac{1}{2}\gamma_{x'y'})$"
        second = r"$(\varepsilon_{y'y'},\,-\frac{1}{2}\gamma_{x'y'})$"

    # Draw the circle.
    sigma = (sxx + syy) / 2 + (sxx - syy) / 2 * \
        np.cos(2 * par2) + txy * np.sin(2 * par2)
    tau = - (sxx - syy) / 2 * np.sin(2 * par2) + txy * np.cos(2 * par2)

    # Draw the circle.
    ax.plot(sigma, -1 ** ccw * tau, color='white', linewidth=2, zorder=1)

    dx = 6
    dy = 3

    # Draw point (sigma_xx, -tau_xy)
    ax.scatter([sxx2], [(-1) ** ccw * txy2], color="yellow", s=50, zorder=5)
    ax.scatter([sxx], [(-1) ** ccw * txy], color="blue", s=50, zorder=5)
    ax.text(sxx2 + 12, (-1) ** ccw * txy2 + 2.9, 'A')

    # Draw point (sigma_yy, tau_xy)
    ax.scatter([syy2], [(-1) ** (ccw+1) * txy2],
               color="yellow", s=50, zorder=5)
    ax.scatter([syy], [(-1) ** (ccw+1) * txy], color="blue", s=50, zorder=5)
    ax.text(syy2 - 10, (-1) ** (ccw+1) * txy2 + 2.9, 'B')

    # Draw line.
    ax.plot([sxx, syy], [(-1) ** ccw * txy, (-1) ** (ccw+1) * txy],
            color='deepskyblue', zorder=1)
    ax.plot([sxx2, syy2], [(-1) ** ccw * txy2,
                           (-1) ** (ccw+1) * txy2], color='red')

    # Draw reference points
    ax.text((sxx + syy) / 2 - dx, dy, 'C')
    ax.scatter((sxx + syy) / 2, 0, color="cyan", s=30, zorder=5)

    ax.text(130.9 + dx/2, dy, 'D')
    ax.text(130.9 + 0.8 * dx, -2.4 * dy, r'$\varepsilon_1$')
    ax.scatter(130.9, 0, color="cyan", s=30, zorder=5)

    ax.text(19.1 + dx/2, dy, 'E')
    ax.text(19.1 - 1.8 * dx, -2.4 * dy, r'$\varepsilon_2$')
    ax.scatter(19.1, 0, color="cyan", s=30, zorder=5)

    # Max shear stress.
    ax.scatter((sxx+syy)/2, 55.9, color='violet', s=30, zorder=7)
    ax.text(-4 * dx, 55, r'$\frac{1}{2}\gamma_{\mathrm{max}}$')
    ax.plot([0, (sxx+syy)/2], [55.9, 55.9], color='white', linestyle='dashed')

    ax.scatter((sxx+syy)/2, -55.9, color='violet', s=30, zorder=7)
    ax.text(-5 * dx, -56, r'$-\frac{1}{2}\gamma_{\mathrm{max}}$')
    ax.plot([0, (sxx+syy)/2], [-55.9, -55.9],
            color='white', linestyle='dashed')

    ax.set_xlim(-30, 270)
    ax.set_ylim(-80, 80)

    # Draw rotating axis.
    ax.text(192, 42, 'O')
    ax.text(-8, -8, 'O')
    ax.text(251, 50, r'$x$')
    ax.text(201, 100, r'$y$')
    ax.arrow(200, 50, 50, 0., fc='deepskyblue', ec='deepskyblue', lw=lw, head_width=hw, head_length=hl,
             overhang=ohg, length_includes_head=True, clip_on=False)
    ax.arrow(200, 50, 0, 50., fc='deepskyblue', ec='deepskyblue', lw=lw, head_width=hw, head_length=hl,
             overhang=ohg, length_includes_head=True, clip_on=False)

    ax.arrow(200, 50, 50*cos(par[i]), 50*sin(par[i]), fc='r', ec='r', lw=lw, head_width=hw, head_length=hl,
             overhang=ohg, length_includes_head=True, clip_on=False)
    ax.arrow(200, 50, -50*sin(par[i]), 50*cos(par[i]), fc='r', ec='r', lw=lw, head_width=hw, head_length=hl,
             overhang=ohg, length_includes_head=True, clip_on=False)

    ax.text(200 + 50*cos(par[i]) + 1, 50 + 50*sin(par[i]), r"$x'$")
    ax.text(200 - 50*sin(par[i]) + 1, 50 + 50*cos(par[i]), r"$y'$")

    # Draw labels.
    ax.text(200, 0, r'A =' + first)
    ax.text(200, -20, r'B =' + second)
    ax.text(
        200, -40, r'C = $\big(\frac{\varepsilon_1+\varepsilon_2}{2}, 0\big)$')
    ax.text(200, -60, r'D = $\big(\varepsilon_1, 0\big)$')
    ax.text(200, -80, r'E = $\big(\varepsilon_2, 0\big)$')
    ax.text(200, -100, r'$\gamma_{\mathrm{max}}=\varepsilon_1-\varepsilon_2$')

    # Draw angle.
    if ccw == 1:
        dteta = Wedge(((sxx+syy) / 2, 0), 55.90 / 2, -2 * 31.72, -
                      2 * 31.72 + 2*par[i] * 180 / pi, color='lightgray')
    elif ccw == 0:
        dteta = Wedge(((sxx+syy) / 2, 0), 55.90 / 2, 2 * 31.72 -
                      2*par[i] * 180 / pi, 2 * 31.72, color='lightgray')
    dtetaref = Wedge((200, 50), 25, 0, 180 * par[i] / pi, color='lightgray')
    ax.add_patch(dteta)
    ax.add_patch(dtetaref)

    # Draw theta.
    tp = 31.72 * pi / 180
    ct = (tp - par[i]) + par[i] / 2
    if ccw == 1:
        xt = (sxx + syy) / 2 + 55.9 / 1.5 * cos(2*ct)
        yt = - 55.9 / 1.5 * sin(2*ct)
    else:
        xt = (sxx + syy) / 2 + 55.9 / 1.5 * cos(-2*ct)
        yt = - 55.9 / 1.5 * sin(-2*ct)

    xta = 200 + 30 * cos(par[i] / 2)
    yta = 50 + 30 * sin(par[i] / 2)
    if i == len(par) - 1:
        if ccw == 1:
            ax.text(xt, yt, r'$2\theta_p$')
            ax.text(xta, yta, r'$\theta_p$')
        else:
            ax.text(xt, yt, r'$2\theta_p$')
            ax.text(xta, yta, r'$\theta_p$')
    else:
        ax.text(xt, yt, r'$2\theta$')
        ax.text(xta, yta, r'$\theta$')


# Create the animation.
anim = animation.FuncAnimation(fig, animate, init_func=init,
                               frames=len(par), interval=50, blit=False, repeat=False)

if save:
    Writer = animation.writers['ffmpeg']
    writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=1800)
    anim.save('mohrCircleStrain.mp4', writer=writer)

plt.show()

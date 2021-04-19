# Master in Mechanical Engineering
# Solid Mechanics
#                                   Faculty of Engineering of University of Porto, 2019/2020
# ==========================================================================================
# Problem 42: Torsion of a trianglular prsimatic shaft using Saint-Venant Theory
#             Contour plots
#                                           A. M. Couto Carneiro <amcc@fe.up.pt>, @FEUP/CM2S
# ==========================================================================================
# Import modules
# --------------
# Array calcuations
import numpy as np
# Plotting toolbox
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
from matplotlib import gridspec
# Triangulations
from matplotlib.tri import Triangulation, TriAnalyzer, UniformTriRefiner
# Geometry operations
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon

# Settings
# --------
plt.rc('text', usetex=True)
plt.rc('text.latex', preamble=r'\usepackage{amsmath} \renewcommand{\familydefault}{\sfdefault} \usepackage[cm]{sfmath}')
plt.rc('xtick', direction="in")
plt.rc('ytick', direction="in")
plt.rcParams.update({'font.size': 20})

# Functions
# ---------
def isInsideTri(x, y):
    """Checks if a point is inside the triangle.
    """
    point = Point(x, y)
    return triangle.contains(point)

def Phi(x, y):
    """Saint-Venant torsion function.
    """
    Phi = 1 / 2 * (x - np.sqrt(3) * y - 2 / 3 * h) * (x + np.sqrt(3) * y - 2 / 3 * h) * (x + h / 3)
    return Phi

def shearStress(x, y):
    """Shear stress componenents.
    """
    tauxz = 1 / h * (3 * y * x + h * y)
    tauyz = 1 / (2 * h) * (3 * x ** 2 - 2 * h * x - 3 * y ** 2)
    return tauxz, tauyz
# ==========================================================================================
# Number of random points to be generated
N = 1000
# Number of equally spaced points over the edges
nEdge = 30
# Height of the triangle
h = 1
# Triangle definition
triangle = Polygon([(-h / 3, -h / np.sqrt(3)), (2 * h / 3, 0), (- h / 3, h / np.sqrt(3))])
# Minimum radius for an inscribed circle on a triangle
minRadius = (h / nEdge) / 10
# Number of subdivisions for triangulation refinement
subdiv = 3

# Generate random points inside a rectangle containing the triangle
randomGen = np.random.RandomState(seed=1)
xTest = randomGen.uniform(- h / 3, 2 * h / 3, size=N)
yTest = randomGen.uniform(-h / np.sqrt(3), h / np.sqrt(3), size=N)

# Take only the points inside the triangle
insidePoints = np.vectorize(isInsideTri)(xTest, yTest)
x = xTest[insidePoints]
y = yTest[insidePoints]

# Inlcude the points on the boundaries in the triangulation
x = np.append(x, np.linspace(- h / 3, 2 * h / 3, num = nEdge, endpoint=False))
x = np.append(x, np.linspace(2 * h / 3, - h / 3, num = nEdge, endpoint=False))
x = np.append(x, np.linspace(- h / 3, - h / 3, num = nEdge, endpoint=False))

y = np.append(y, np.linspace(- h / np.sqrt(3), 0, num = nEdge, endpoint=False))
y = np.append(y, np.linspace(0, h / np.sqrt(3), num = nEdge, endpoint=False))
y = np.append(y, np.linspace(h / np.sqrt(3), -h / np.sqrt(3), num = nEdge, endpoint=False))

# Triangulate the scattered points
triMesh = Triangulation(x, y)

# Remove bad triangles from plot
mask = TriAnalyzer(triMesh).get_flat_tri_mask(minRadius)
# triMesh.set_mask(mask)

# DEBUG: Bad triangles in mesh
badTriMesh = Triangulation(x, y)
badTriMesh.set_mask(~mask)

# Compute the height values
phiVal = Phi(triMesh.x, triMesh.y)
tauxz, tauyz = shearStress(triMesh.x, triMesh.y)

# Refine the trianglulation
refiner = UniformTriRefiner(triMesh)
triMeshRef = refiner.refine_field(phiVal, subdiv=subdiv)
triMeshRef, tauxzRef = refiner.refine_field(tauxz, subdiv=subdiv)
triMeshRef, tauyzRef = refiner.refine_field(tauyz, subdiv=subdiv)
phiValRef = Phi(triMeshRef.x, triMeshRef.y)
tauRef = np.sqrt(tauxzRef ** 2 + tauyzRef ** 2)

# Plot the triangulations
# =======================
fig1, ax1 = plt.subplots(1, 1, num="Triangulations", figsize=(8.5, 8), constrained_layout=True)
ax1.set_aspect("equal")
ax1.triplot(triMeshRef, color='gray', linewidth=0.5)
ax1.triplot(triMesh, color='blue', linewidth=1.0)
ax1.triplot(badTriMesh, color='red', linewidth=2.0)

refLine = ax1.plot([0, 0], [0, 0], color='gray', linewidth=0.5, label="Refined")
delLine = ax1.plot([0, 0], [0, 0], color='blue', linewidth=1.0, label="Delaunay")
badLine = ax1.plot([0, 0], [0, 0], color='red', linewidth=2.0, label="Bad Triangles")
ax1.legend(loc="upper right", frameon=False)
ax1.axis("off")

# Plot the potential function
# ===========================
fig2 = plt.figure(num="Saint-Venant Function", figsize=(8.5, 8), constrained_layout=False)
spec = fig2.add_gridspec(ncols=2, nrows=1, bottom=0.02, top=0.98, left=0.15, right=1.0, width_ratios=[0.02, 0.9], height_ratios=[1.0], hspace=0.0)
ax2 = fig2.add_subplot(spec[1])
ax2.set_aspect("equal")
minval = np.min(phiValRef)
maxval = np.max(phiValRef)
levels = MaxNLocator(nbins=10).tick_values(minval, maxval)
cmap = plt.get_cmap('viridis')
norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
cplot = ax2.tricontourf(triMeshRef, phiValRef, alpha=0.85, levels=levels, cmap=cmap, norm=norm)
ax2.tricontour(triMeshRef, phiValRef, levels=levels, cmap=cmap, linewidths=[2.0], norm=norm)
ax2.axis("off")
cax2 = fig2.add_subplot(spec[0])
cbar = plt.colorbar(cplot, cax=cax2, orientation="vertical")
cax2.yaxis.set_label_position("left")
cax2.yaxis.tick_left()
cax2.set_ylabel(r"$\dfrac{\Phi}{G\theta}$", rotation=0, labelpad=20)
ax2.plot([-h/3, 2*h/3, -h/3, -h/3], [-h/np.sqrt(3), 0, h/np.sqrt(3), -h/np.sqrt(3)], color='black', linewidth=2.0)
ax2.set_xlim([-h/3 - h/100, 2*h/3 + h / 100])
ax2.set_ylim([-h/np.sqrt(3) - h/100, h/np.sqrt(3) + h / 100])
fig2.savefig("saintVenantFunction.pdf", format="pdf", bbox_inches="tight")

# Plot the xz shear stress
# =========================
fig3 = plt.figure(num="Shear stress xz", figsize=(8.5, 8), constrained_layout=False)
spec = fig3.add_gridspec(ncols=2, nrows=1, bottom=0.02, top=0.98, left=0.15, right=1.0, width_ratios=[0.02, 0.9], height_ratios=[1.0], hspace=0.0)
ax3 = fig3.add_subplot(spec[1])
ax3.set_aspect("equal")
minval = np.min(tauxzRef)
maxval = np.max(tauxzRef)
levels = MaxNLocator(nbins=15).tick_values(-0.5, +0.5)
cmap = plt.get_cmap('magma')
norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
cplot = ax3.tricontourf(triMeshRef, tauxzRef, alpha=0.85, levels=levels, cmap=cmap, norm=norm)
ax3.tricontour(triMeshRef, tauxzRef, levels=levels, cmap=cmap, linewidths=[2.0], norm=norm)
ax3.axis("off")
cax3 = fig3.add_subplot(spec[0])
cbar = plt.colorbar(cplot, cax=cax3, orientation="vertical")
cax3.yaxis.set_label_position("left")
cax3.yaxis.tick_left()
cax3.set_ylabel(r"$\dfrac{\tau_{xz}}{G\theta h}$", rotation=0, labelpad=20)
ax3.plot([-h/3, 2*h/3, -h/3, -h/3], [-h/np.sqrt(3), 0, h/np.sqrt(3), -h/np.sqrt(3)], color='black', linewidth=2.0)
ax3.set_xlim([-h/3 - h/100, 2*h/3 + h / 100])
ax3.set_ylim([-h/np.sqrt(3) - h/100, h/np.sqrt(3) + h / 100])
fig3.savefig("stressxz.pdf", format="pdf", bbox_inches="tight")

# Plot the yz shear stress
# =========================f
fig4 = plt.figure(num="Shear stress yz", figsize=(8.5, 8), constrained_layout=False)
spec = fig4.add_gridspec(ncols=2, nrows=1, bottom=0.02, top=0.98, left=0.15, right=1.0, width_ratios=[0.02, 0.9], height_ratios=[1.0], hspace=0.0)
ax4 = fig4.add_subplot(spec[1])
ax4.set_aspect("equal")
minval = np.min(tauyzRef)
maxval = np.max(tauyzRef)
levels = MaxNLocator(nbins=15).tick_values(-0.32, +0.5)
cmap = plt.get_cmap('cividis')
norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
cplot = ax4.tricontourf(triMeshRef, tauyzRef, alpha=0.85, levels=levels, cmap=cmap, norm=norm)
ax4.tricontour(triMeshRef, tauyzRef, levels=levels, cmap=cmap, linewidths=[2.0], norm=norm)
ax4.axis("off")
cax4 = fig4.add_subplot(spec[0])
cbar = plt.colorbar(cplot, cax=cax4, orientation="vertical")
cax4.yaxis.set_label_position("left")
cax4.yaxis.tick_left()
cax4.set_ylabel(r"$\dfrac{\tau_{yz}}{G\theta h}$", rotation=0, labelpad=20)
ax4.plot([-h/3, 2*h/3, -h/3, -h/3], [-h/np.sqrt(3), 0, h/np.sqrt(3), -h/np.sqrt(3)], color='black', linewidth=2.0)
ax4.set_xlim([-h/3 - h/100, 2*h/3 + h / 100])
ax4.set_ylim([-h/np.sqrt(3) - h/100, h/np.sqrt(3) + h / 100])
fig4.savefig("stressyz.pdf", format="pdf", bbox_inches="tight")

# Plot the resulting shear stress
# ===============================
fig5 = plt.figure(num="Resultant shear stress", figsize=(8.5, 8), constrained_layout=False)
spec = fig5.add_gridspec(ncols=2, nrows=1, bottom=0.02, top=0.98, left=0.15, right=1.0, width_ratios=[0.02, 0.9], height_ratios=[1.0], hspace=0.0)
ax5 = fig5.add_subplot(spec[1])
ax5.set_aspect("equal")
minval = np.min(tauRef)
maxval = np.max(tauRef)
levels = MaxNLocator(nbins=10).tick_values(0, 0.5)
cmap = plt.get_cmap('plasma')
norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
cplot = ax5.tricontourf(triMeshRef, tauRef, alpha=0.85, levels=levels, cmap=cmap, norm=norm)
ax5.tricontour(triMeshRef, tauRef, levels=levels, cmap=cmap, linewidths=[2.0], norm=norm)
ax5.axis("off")
cax5 = fig5.add_subplot(spec[0])
cbar = plt.colorbar(cplot, cax=cax5, orientation="vertical")
cax5.yaxis.set_label_position("left")
cax5.yaxis.tick_left()
cax5.set_ylabel(r"$\dfrac{\tau}{G \theta h}$", rotation=0, labelpad=20)
ax5.plot([-h/3, 2*h/3, -h/3, -h/3], [-h/np.sqrt(3), 0, h/np.sqrt(3), -h/np.sqrt(3)], color='black', linewidth=2.0)
ax5.set_xlim([-h/3 - h/100, 2*h/3 + h / 100])
ax5.set_ylim([-h/np.sqrt(3) - h/100, h/np.sqrt(3) + h / 100])
fig5.savefig("resultantShear.pdf", format="pdf", bbox_inches="tight")

plt.show()

import gmsh
import sys

# Initialize Gmsh API
gmsh.initialize()

# Suppress verbose messages
gmsh.option.setNumber("General.Terminal", 1)

# Add a new model named '3D Box'
gmsh.model.add("3D Box")

# Define the box dimensions
length = 1.0
width = 1.0
height = 1.0
lc = 12

# Create a 3D box geometry (Gmsh entities: Point, Line, Surface, Volume)
# Define the 8 corner points of the box
p1 = gmsh.model.geo.addPoint(0, 0, 0, lc)
p2 = gmsh.model.geo.addPoint(length, 0, 0, lc)
p3 = gmsh.model.geo.addPoint(length, width, 0, lc)
p4 = gmsh.model.geo.addPoint(0, width, 0, lc)
p5 = gmsh.model.geo.addPoint(0, 0, height, lc)
p6 = gmsh.model.geo.addPoint(length, 0, height, lc)
p7 = gmsh.model.geo.addPoint(length, width, height, lc)
p8 = gmsh.model.geo.addPoint(0, width, height, lc)

# Create the 12 edges (lines) connecting the points
l1 = gmsh.model.geo.addLine(p1, p2)
l2 = gmsh.model.geo.addLine(p2, p3)
l3 = gmsh.model.geo.addLine(p3, p4)
l4 = gmsh.model.geo.addLine(p4, p1)
l5 = gmsh.model.geo.addLine(p5, p6)
l6 = gmsh.model.geo.addLine(p6, p7)
l7 = gmsh.model.geo.addLine(p7, p8)
l8 = gmsh.model.geo.addLine(p8, p5)
l9 = gmsh.model.geo.addLine(p1, p5)
l10 = gmsh.model.geo.addLine(p2, p6)
l11 = gmsh.model.geo.addLine(p3, p7)
l12 = gmsh.model.geo.addLine(p4, p8)

# Create the six faces (loops of four lines)
cl1 = gmsh.model.geo.addCurveLoop([l1, l2, l3, l4])
cl2 = gmsh.model.geo.addCurveLoop([l5, l6, l7, l8])
cl3 = gmsh.model.geo.addCurveLoop([l1, l10, -l5, -l9])
cl4 = gmsh.model.geo.addCurveLoop([l2, l11, -l6, -l10])
cl5 = gmsh.model.geo.addCurveLoop([l3, l12, -l7, -l11])
cl6 = gmsh.model.geo.addCurveLoop([l4, l9, -l8, -l12])

f1 = gmsh.model.geo.addPlaneSurface([cl1])
f2 = gmsh.model.geo.addPlaneSurface([cl2])
f3 = gmsh.model.geo.addPlaneSurface([cl3])
f4 = gmsh.model.geo.addPlaneSurface([cl4])
f5 = gmsh.model.geo.addPlaneSurface([cl5])
f6 = gmsh.model.geo.addPlaneSurface([cl6])

surLoop = gmsh.model.geo.addSurfaceLoop([f1, f2, f3, f4, f5, f6])
# Define the volume of the box using the 6 faces
volume = gmsh.model.geo.addVolume([surLoop])

# Synchronize the geometry
gmsh.model.geo.synchronize()

# Specify the mesh element size
element_size = 0.1  # Adjust the size to control mesh resolution

# Generate 3D tetrahedral mesh for the volume
gmsh.model.mesh.generate(3)

# Save the mesh to a file (can save in .msh format or others like .vtk)
gmsh.write("box_tets_3D.m")

# Launch Gmsh's graphical interface to view the geometry and mesh
if '-nopopup' not in sys.argv:
    gmsh.fltk.run()

# Finalize the Gmsh session
gmsh.finalize()

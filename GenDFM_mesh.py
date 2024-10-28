import gmsh
import sys

def addVolume(lc):
    length = 1.0
    width = 1.0
    height = 1.0

    # Create a 3D box geometry (Gmsh entities: Point, Line, Surface, Volume)
    # Define the 8 corner points of the box
    p1 = gmsh.model.occ.addPoint(0, 0, 0, lc)
    p2 = gmsh.model.occ.addPoint(length, 0, 0, lc)
    p3 = gmsh.model.occ.addPoint(length, width, 0, lc)
    p4 = gmsh.model.occ.addPoint(0, width, 0, lc)
    p5 = gmsh.model.occ.addPoint(0, 0, height, lc)
    p6 = gmsh.model.occ.addPoint(length, 0, height, lc)
    p7 = gmsh.model.occ.addPoint(length, width, height, lc)
    p8 = gmsh.model.occ.addPoint(0, width, height, lc)

    # Create the 12 edges (lines) connecting the points
    l1 = gmsh.model.occ.addLine(p1, p2)
    l2 = gmsh.model.occ.addLine(p2, p3)
    l3 = gmsh.model.occ.addLine(p3, p4)
    l4 = gmsh.model.occ.addLine(p4, p1)
    l5 = gmsh.model.occ.addLine(p5, p6)
    l6 = gmsh.model.occ.addLine(p6, p7)
    l7 = gmsh.model.occ.addLine(p7, p8)
    l8 = gmsh.model.occ.addLine(p8, p5)
    l9 = gmsh.model.occ.addLine(p1, p5)
    l10 = gmsh.model.occ.addLine(p2, p6)
    l11 = gmsh.model.occ.addLine(p3, p7)
    l12 = gmsh.model.occ.addLine(p4, p8)

    # Create the six faces (loops of four lines)
    cl1 = gmsh.model.occ.addCurveLoop([l1, l2, l3, l4])
    cl2 = gmsh.model.occ.addCurveLoop([l5, l6, l7, l8])
    cl3 = gmsh.model.occ.addCurveLoop([l1, l10, -l5, -l9])
    cl4 = gmsh.model.occ.addCurveLoop([l2, l11, -l6, -l10])
    cl5 = gmsh.model.occ.addCurveLoop([l3, l12, -l7, -l11])
    cl6 = gmsh.model.occ.addCurveLoop([l4, l9, -l8, -l12])

    f1 = gmsh.model.occ.addPlaneSurface([cl1])
    f2 = gmsh.model.occ.addPlaneSurface([cl2])
    f3 = gmsh.model.occ.addPlaneSurface([cl3])
    f4 = gmsh.model.occ.addPlaneSurface([cl4])
    f5 = gmsh.model.occ.addPlaneSurface([cl5])
    f6 = gmsh.model.occ.addPlaneSurface([cl6])

    surLoop = gmsh.model.occ.addSurfaceLoop([f1, f2, f3, f4, f5, f6])
    # Define the volume of the box using the 6 faces
    volume = gmsh.model.occ.addVolume([surLoop])
    return volume


# Initialize Gmsh API
gmsh.initialize()

# Suppress verbose messages
gmsh.option.setNumber("General.Terminal", 1)

gmsh.model.add("DFM")


lc = 0.15

fc_p1 = gmsh.model.occ.addPoint(0.5, 0, 0, lc)
fc_p2 = gmsh.model.occ.addPoint(0.5, 0, 1.5, lc)
fc_p3 = gmsh.model.occ.addPoint(0.5, 1.5, 1.5, lc)
fc_p4 = gmsh.model.occ.addPoint(0.5, 1.5, 0, lc)

fc_l1 = gmsh.model.occ.addLine(fc_p1, fc_p2)
fc_l2 = gmsh.model.occ.addLine(fc_p2, fc_p3)
fc_l3 = gmsh.model.occ.addLine(fc_p3, fc_p4)
fc_l4 = gmsh.model.occ.addLine(fc_p4, fc_p1)

fc_p5 = gmsh.model.occ.addPoint(1, 1, 0.5, lc)
fc_p6 = gmsh.model.occ.addPoint(1, 0, 0.5, lc)
fc_p7 = gmsh.model.occ.addPoint(0, 0, 0.5, lc)
fc_p8 = gmsh.model.occ.addPoint(0, 1, 0.5, lc)

fc_l5 = gmsh.model.occ.addLine(fc_p5, fc_p6)
fc_l6 = gmsh.model.occ.addLine(fc_p6, fc_p7)
fc_l7 = gmsh.model.occ.addLine(fc_p7, fc_p8)
fc_l8 = gmsh.model.occ.addLine(fc_p8, fc_p5)

fc_cl1 = gmsh.model.occ.addCurveLoop([fc_l1, fc_l2, fc_l3, fc_l4])
fc_cl2 = gmsh.model.occ.addCurveLoop([fc_l5, fc_l6,fc_l7, fc_l8])

fc_f1 = gmsh.model.occ.addPlaneSurface([fc_cl1])
fc_f2 = gmsh.model.occ.addPlaneSurface([fc_cl2])
gmsh.model.occ.synchronize()


gmsh.model.occ.synchronize()

volumeN1 = addVolume(lc)

outDimTags, outDimTagMap = gmsh.model.occ.intersect([(3, volumeN1)], [(2, fc_f1), (2, fc_f2)])
#gmsh.model.occ.intersect([(3, volume)], [(2, fc_f2)])

gmsh.model.occ.synchronize()

mesh_size = lc

gmsh.option.setNumber("Mesh.MeshSizeMin", lc)
gmsh.option.setNumber("Mesh.MeshSizeMax", mesh_size)
gmsh.option.setNumber("Mesh.Algorithm", 1)

gmsh.model.mesh.generate(2)
if "-nopopup" not in sys.argv:
    gmsh.fltk.run()

volumeN2 = addVolume(lc)
gmsh.model.occ.synchronize()

gmsh.model.occ.fragment([(3, volumeN2)], outDimTags)
gmsh.model.occ.synchronize()

gmsh.model.mesh.generate(3)
if "-nopopup" not in sys.argv:
    gmsh.fltk.run()

gmsh.write("DFM_mesh.m")

gmsh.finalize()

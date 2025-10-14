import gmsh
import sys

gmsh.initialize()
gmsh.model.add("SleeperExtrusion")

sleeper_profile_points = [
    (-0.756813, 0.172471),
    (-0.756813, 0.061628),
    (-0.312048, 0.061628),
    (-0.312048, 0.031628),
    (-0.132048, 0.031628),
    (-0.132048, 0.061628),
    (1.189789, 0.061628),
    (1.189789, 0.031628),
    (1.369789, 0.031628),
    (1.369789, 0.061628),
    (1.843187, 0.061628),
    (1.843187, 0.161628),
    (1.543187, 0.261628),
    (1.369789, 0.261628),
    (1.369789, 0.281628),
    (1.189789, 0.281628),
    (1.189789, 0.261628),
    (1.010661, 0.261628),
    (0.710661, 0.161628),
    (0.310661, 0.161628),
    (0.010661, 0.261628),
    (-0.132048, 0.261628),
    (-0.132048, 0.281628),
    (-0.312048, 0.281628),
    (-0.312048, 0.261628),
    (-0.489339, 0.261628),
]

sleeper_y_coords = [y for y, _ in sleeper_profile_points]
sleeper_z_coords = [z for _, z in sleeper_profile_points]
sleeper_y_center = (min(sleeper_y_coords) + max(sleeper_y_coords)) / 2
sleeper_z_center = (min(sleeper_z_coords) + max(sleeper_z_coords)) / 2

sleeper_profile_points = [(y - sleeper_y_center, z - sleeper_z_center) for y, z in sleeper_profile_points]

sleeper_profile_pts = []
for i, (y, z) in enumerate(sleeper_profile_points):
    sleeper_profile_pts.append(gmsh.model.occ.addPoint(0, y, z))

sleeper_profile_lines = []
for i in range(len(sleeper_profile_pts)):
    start_idx = i
    end_idx = (i + 1) % len(sleeper_profile_pts)
    sleeper_profile_lines.append(gmsh.model.occ.addLine(sleeper_profile_pts[start_idx], sleeper_profile_pts[end_idx]))

sleeper_cl = gmsh.model.occ.addCurveLoop(sleeper_profile_lines)
sleeper_profile_face = gmsh.model.occ.addPlaneSurface([sleeper_cl])

extrusion_length = 0.3
start_point = gmsh.model.occ.addPoint(0, 0, 0)
end_point = gmsh.model.occ.addPoint(extrusion_length, 0, 0)
extrusion_path = gmsh.model.occ.addLine(start_point, end_point)
gmsh.model.occ.addWire([extrusion_path], 3)
sleeper_sweep = gmsh.model.occ.addPipe([(2, sleeper_profile_face)], 3)

gmsh.model.occ.synchronize()
gmsh.write("sleeper_sweep.step")
gmsh.finalize()

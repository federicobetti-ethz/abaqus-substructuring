import gmsh
import sys


def generate_sleeper_sweep(
    sleeper_length: float,
    sleeper_profile_points: list,
    step_file: str = "sleeper_sweep.step",
    model_name: str = "SleeperExtrusion",
):
    gmsh.initialize()
    gmsh.model.add(model_name)

    sleeper_y_coords = [y for y, _ in sleeper_profile_points]
    sleeper_z_coords = [z for _, z in sleeper_profile_points]
    sleeper_y_center = (min(sleeper_y_coords) + max(sleeper_y_coords)) / 2
    sleeper_z_center = (min(sleeper_z_coords) + max(sleeper_z_coords)) / 2

    sleeper_profile_points = [
        (y - sleeper_y_center, z - sleeper_z_center) for y, z in sleeper_profile_points
    ]

    sleeper_profile_pts = []
    for i, (y, z) in enumerate(sleeper_profile_points):
        sleeper_profile_pts.append(gmsh.model.occ.addPoint(0, -y, -z))

    sleeper_profile_lines = []
    for i in range(len(sleeper_profile_pts)):
        start_idx = i
        end_idx = (i + 1) % len(sleeper_profile_pts)
        sleeper_profile_lines.append(
            gmsh.model.occ.addLine(
                sleeper_profile_pts[start_idx], sleeper_profile_pts[end_idx]
            )
        )

    sleeper_cl = gmsh.model.occ.addCurveLoop(sleeper_profile_lines)
    sleeper_profile_face = gmsh.model.occ.addPlaneSurface([sleeper_cl])

    start_point = gmsh.model.occ.addPoint(0, 0, 0)
    end_point = gmsh.model.occ.addPoint(sleeper_length, 0, 0)
    extrusion_path = gmsh.model.occ.addLine(start_point, end_point)
    gmsh.model.occ.addWire([extrusion_path], 3)
    sleeper_sweep = gmsh.model.occ.addPipe([(2, sleeper_profile_face)], 3)

    gmsh.model.occ.synchronize()
    gmsh.write(step_file)
    gmsh.finalize()

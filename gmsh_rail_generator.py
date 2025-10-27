import gmsh

def generate_rail_sweep(txt_file: str, rail_profile_points: list):
    step_file = "rail_sweep.step"
    
    gmsh.initialize()
    gmsh.model.add("RailSweep")

    y_coords = [y for y, _ in rail_profile_points]
    z_coords = [z for _, z in rail_profile_points]
    y_center = (min(y_coords) + max(y_coords)) / 2
    z_center = (min(z_coords) + max(z_coords)) / 2

    rail_profile_points = [(y - y_center, z - z_center) for y, z in rail_profile_points]

    profile_pts = []
    for i, (y, z) in enumerate(rail_profile_points):
        profile_pts.append(gmsh.model.occ.addPoint(0, -y, -z))

    profile_lines = []
    for i in range(len(profile_pts)):
        start_idx = i
        end_idx = (i + 1) % len(profile_pts)
        profile_lines.append(gmsh.model.occ.addLine(profile_pts[start_idx], profile_pts[end_idx]))

    cl = gmsh.model.occ.addCurveLoop(profile_lines)
    profile_face = gmsh.model.occ.addPlaneSurface([cl])

    ss, xs, ys, zs = [], [], [], []
    with open(txt_file, "r") as f:
        lines = f.readlines()
        for idx, line in enumerate(lines):
            if line.startswith('"Superelevation u(s)"'):
                superelevation_line = idx
            elif line.startswith('"x(s)"'):
                xline = idx
            elif line.startswith('"y(s)"'):
                yline = idx
            elif line.startswith('"z(s)"'):
                zline = idx
        
        for line in lines[superelevation_line+3:xline]:
            ss.append(float(line.split(",")[0]))
        for line in lines[xline+3:yline]:
            xs.append(float(line.split(",")[1]))
        for line in lines[yline+3:zline]:
            ys.append(float(line.split(",")[1]))
        for line in lines[zline+3:]:
            zs.append(float(line.split(",")[1]))

    spline_pts = [gmsh.model.occ.addPoint(xs[i], ys[i], zs[i]) for i in range(len(xs))]
    spline_pts = spline_pts[::10]
    spline = gmsh.model.occ.addBSpline(spline_pts)
    gmsh.model.occ.addWire([spline], 3)
    sweep = gmsh.model.occ.addPipe([(2, profile_face)], 3)
    gmsh.model.occ.synchronize()
    gmsh.write(step_file)
    gmsh.finalize()

    return step_file
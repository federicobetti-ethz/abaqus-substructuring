from abaqus import *
from abaqusConstants import *
from part import *
from material import *
from section import *
from assembly import *
from step import *
from interaction import *
from load import *
from mesh import *
from optimization import *
from job import *
from sketch import *
from visualization import *
from connectorBehavior import *

import numpy as np
from scipy.interpolate import CubicSpline, interp1d

EPS = 1e-2

rail_profile_points = [
    (1.316689, 423.774E-03),
    (1.307883, 418.203E-03),
    (1.297689, 414.324E-03),
    (1.293897, 410.876E-03),
    (1.291245, 403.458E-03),
    (1.289639, 395.455E-03),
    (1.288447, 385.816E-03),
    (1.288039, 375.928E-03),
    (1.288039, 365.147E-03),
    (1.288039, 354.478E-03),
    (1.288039, 343.928E-03),
    (1.28855, 332.865E-03),
    (1.290056, 322.022E-03),
    (1.292114, 313.325E-03),
    (1.29636, 308.967E-03),
    (1.304345, 306.064E-03),
    (1.314243, 302.464E-03),
    (1.324179, 298.851E-03),
    (1.334994, 296.545E-03),
    (1.351075, 295.394E-03),
    (1.354789, 291.404E-03),
    (1.354789, 285.628E-03),
    (1.340813, 283.628E-03),
    (1.331635, 283.628E-03),
    (1.320693, 283.628E-03),
    (1.310722, 283.628E-03),
    (1.299821, 283.628E-03),
    (1.289311, 283.628E-03),
    (1.279789, 283.628E-03),
    (1.270053, 283.628E-03),
    (1.259174, 283.628E-03),
    (1.247651, 283.628E-03),
    (1.236811, 283.628E-03),
    (1.227038, 283.628E-03),
    (1.217987, 283.628E-03),
    (1.206789, 283.628E-03),
    (1.204789, 291.404E-03),
    (1.208504, 295.394E-03),
    (1.217572, 296.043E-03),
    (1.224584, 296.545E-03),
    (1.235399, 298.851E-03),
    (1.244588, 302.192E-03),
    (1.254797, 305.905E-03),
    (1.263218, 308.967E-03),
    (1.267464, 313.325E-03),
    (1.269523, 322.022E-03),
    (1.271028, 332.863E-03),
    (1.271539, 343.928E-03),
    (1.271539, 354.478E-03),
    (1.271539, 365.147E-03),
    (1.271539, 375.928E-03),
    (1.271131, 385.814E-03),
    (1.269941, 395.449E-03),
    (1.268339, 403.458E-03),
    (1.265682, 410.873E-03),
    (1.26189, 414.324E-03),
    (1.25177, 418.175E-03),
    (1.244818, 420.82E-03),
    (1.243191, 429.775E-03),
    (1.243773, 441.367E-03),
    (1.246324, 448.472E-03),
    (1.253665, 453.342E-03),
    (1.261358, 454.77E-03),
    (1.269539, 455.453E-03),
    (1.279789, 455.628E-03),
    (1.290039, 455.453E-03),
    (1.29804, 454.795E-03),
    (1.305914, 453.342E-03),
    (1.312847, 448.992E-03),
    (1.315806, 441.367E-03),
    (1.316355, 430.434E-03),
]

y_coords = [y for y, _ in rail_profile_points]
z_coords = [z for _, z in rail_profile_points]
y_center = (min(y_coords) + max(y_coords)) / 2
z_center = (min(z_coords) + max(z_coords)) / 2

rail_profile_points = [(y - y_center, z - z_center) for y, z in rail_profile_points]

def find_top_point_for_wheel_contact(rail_profile_points):
    zs = [z for _, z in rail_profile_points]
    top_z_index = zs.index(max(zs))
    top_point = (0, rail_profile_points[top_z_index][0], rail_profile_points[top_z_index][1])
    return top_point

def find_lower_points_for_coupling(rail_profile_points):
    zs = [z for _, z in rail_profile_points]
    lower_points = []
    for i, (y, z) in enumerate(rail_profile_points):
        if np.isclose(z, min(zs)):
            lower_points.append((0, y, z))
    return lower_points

def find_edge_along_spline_given_initial_point(rail_part, initial_point):
    edge_to_find = None
    edges = rail_part.edges
    vertices = rail_part.vertices
    for edge in edges:
        vertices_indices = edge.getVertices()
        if np.allclose(vertices[vertices_indices[0]].pointOn[0], initial_point) and vertices[vertices_indices[1]].pointOn[0][0] > 0.0:
            edge_to_find = edge.index
            break
        elif np.allclose(vertices[vertices_indices[1]].pointOn[0], initial_point) and vertices[vertices_indices[0]].pointOn[0][0] > 0.0:
            edge_to_find = edge.index
            break
        else:
            continue
    return edge_to_find

def get_spline_coordinates(file_path):
    ss, xs, ys, zs = [], [], [], []

    with open(file_path, "r") as f:
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

        for line in lines[xline+3:yline]:
            ss.append(float(line.split(",")[0]))
            xs.append(float(line.split(",")[1]))
        for line in lines[yline+3:zline]:
            ys.append(float(line.split(",")[1]))
        for line in lines[zline+3:superelevation_line]:
            zs.append(float(line.split(",")[1]))

    return np.array(ss), np.array(xs), np.array(ys), np.array(zs)

txt_file = "C:/Users/bettif/Documents/Github/estra/out/results_0/models/Track.output/Track-Trk_Track.txt"
ss, xs, ys, zs = get_spline_coordinates(txt_file)

xs_spline_interp = interp1d(ss, xs)
ys_spline_interp = interp1d(ss, ys)
zs_spline_interp = interp1d(ss, zs)

mdb = openMdb(pathName='C:/Users/bettif/Documents/GitHub/estra/SampleSplineTrack.cae')

model = mdb.Model(name='SplineBasedTrack', modelType=STANDARD_EXPLICIT)
mdb.openStep('C:/Users/bettif/Documents/GitHub/estra/rail_sweep.step', scaleFromFile=OFF)
model.PartFromGeometryFile(combine=False, dimensionality=
    THREE_D, geometryFile=mdb.acis, name='Rail', type=DEFORMABLE_BODY)

rail_part = model.parts['Rail']
rail_part.seedPart(size=0.05, deviationFactor=0.1, minSizeFactor=0.2)
rail_part.generateMesh()

top_point = find_top_point_for_wheel_contact(rail_profile_points)
lower_points = find_lower_points_for_coupling(rail_profile_points)

edges = rail_part.edges

top_edge = find_edge_along_spline_given_initial_point(rail_part, top_point)
edge_set = edges[top_edge:top_edge+1]
rail_part.Set(name="TopRailEdge", edges=edge_set)
nodes = rail_part.sets["TopRailEdge"].nodes
rail_part.Set(name="TopRailNodes", nodes=nodes)
del rail_part.sets["TopRailEdge"]

lower_nodes = []
for idx, lower_point in enumerate(lower_points):
    lower_edge = find_edge_along_spline_given_initial_point(rail_part, lower_point)
    edge_set = edges[lower_edge:lower_edge+1]
    rail_part.Set(name=f"LowerRailEdge{idx}", edges=edge_set)
    nodes = rail_part.sets[f"LowerRailEdge{idx}"].nodes
    lower_nodes.extend(nodes)
    del rail_part.sets[f"LowerRailEdge{idx}"]

sleeper_spacing = 0.7
sleeper_length = 0.29
rail_burn_in = 0.205

j = 1
while j * sleeper_spacing - rail_burn_in <= xs[-1]:
    labels_nodes_to_keep = []
    xs_start = j * sleeper_spacing - sleeper_length - rail_burn_in
    xs_end = j * sleeper_spacing + sleeper_length - rail_burn_in
    for node in lower_nodes:
        if node.coordinates[0] >= xs_start and node.coordinates[0] <= xs_end:
            labels_nodes_to_keep.append(node.label)
        
    nodes_to_keep = rail_part.nodes.sequenceFromLabels(labels=labels_nodes_to_keep)
    rail_part.Set(name=f"Sleeper{j}", nodes=nodes_to_keep)
    j += 1

mdb.openStep('C:/Users/bettif/Documents/GitHub/estra/sleeper_sweep.step', scaleFromFile=OFF)
model.PartFromGeometryFile(combine=False, dimensionality=
    THREE_D, geometryFile=mdb.acis, name='Sleeper', type=DEFORMABLE_BODY)
sleeper_part = model.parts["Sleeper"]
sleeper_part.seedPart(size=0.1, deviationFactor=0.1, minSizeFactor=0.2)
sleeper_part.generateMesh()

assembly = model.rootAssembly

for i in range(j-1):
    sleeper_instance = assembly.Instance(
        name=f"Sleeper-{i+1}",
        part=sleeper_part,
        dependent=ON,
    )

    s = i * sleeper_spacing + EPS
    x, y, z = xs_spline_interp(s), ys_spline_interp(s), zs_spline_interp(s)

    deriv_y = (ys_spline_interp(s+EPS) - ys_spline_interp(s-EPS)) / (2*EPS)
    deriv_x = (xs_spline_interp(s+EPS) - xs_spline_interp(s-EPS)) / (2*EPS)

    tangent = np.array([deriv_x, deriv_y])
    tangent = tangent / np.linalg.norm(tangent)

    x, y, z = float(x), float(y), float(z)
    sleeper_instance.translate(vector=(x, y, z))

    angle_about_y = np.rad2deg(np.arctan2(tangent[1], tangent[0]))
    sleeper_instance.rotateAboutAxis(
        angle=angle_about_y,
        axisDirection=(0, 0, 1),
        axisPoint=(x, y, z),
    )

rail_left_instance = assembly.Instance(
    name="Rail-Left",
    part=rail_part,
    dependent=ON,
)
rail_left_instance.translate(vector=(0, -0.77, 0.24))

rail_right_instance = assembly.Instance(
    name="Rail-Right",
    part=rail_part,
    dependent=ON,
)
rail_right_instance.translate(vector=(0, 0.74, 0.24))

mdb.saveAs("SplineBasedTrack.cae")

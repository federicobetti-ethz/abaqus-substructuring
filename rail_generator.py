import FreeCAD as App
import FreeCADGui as Gui
import Part

rail_profile_points = [
    (1.316689, 423.774e-03),
    (1.307883, 418.203e-03),
    (1.297689, 414.324e-03),
    (1.293897, 410.876e-03),
    (1.291245, 403.458e-03),
    (1.289639, 395.455e-03),
    (1.288447, 385.816e-03),
    (1.288039, 375.928e-03),
    (1.288039, 365.147e-03),
    (1.288039, 354.478e-03),
    (1.288039, 343.928e-03),
    (1.28855, 332.865e-03),
    (1.290056, 322.022e-03),
    (1.292114, 313.325e-03),
    (1.29636, 308.967e-03),
    (1.304345, 306.064e-03),
    (1.314243, 302.464e-03),
    (1.324179, 298.851e-03),
    (1.334994, 296.545e-03),
    (1.351075, 295.394e-03),
    (1.354789, 291.404e-03),
    (1.354789, 285.628e-03),
    (1.340813, 283.628e-03),
    (1.331635, 283.628e-03),
    (1.320693, 283.628e-03),
    (1.310722, 283.628e-03),
    (1.299821, 283.628e-03),
    (1.289311, 283.628e-03),
    (1.279789, 283.628e-03),
    (1.270053, 283.628e-03),
    (1.259174, 283.628e-03),
    (1.247651, 283.628e-03),
    (1.236811, 283.628e-03),
    (1.227038, 283.628e-03),
    (1.217987, 283.628e-03),
    (1.206789, 283.628e-03),
    (1.204789, 291.404e-03),
    (1.208504, 295.394e-03),
    (1.217572, 296.043e-03),
    (1.224584, 296.545e-03),
    (1.235399, 298.851e-03),
    (1.244588, 302.192e-03),
    (1.254797, 305.905e-03),
    (1.263218, 308.967e-03),
    (1.267464, 313.325e-03),
    (1.269523, 322.022e-03),
    (1.271028, 332.863e-03),
    (1.271539, 343.928e-03),
    (1.271539, 354.478e-03),
    (1.271539, 365.147e-03),
    (1.271539, 375.928e-03),
    (1.271131, 385.814e-03),
    (1.269941, 395.449e-03),
    (1.268339, 403.458e-03),
    (1.265682, 410.873e-03),
    (1.26189, 414.324e-03),
    (1.25177, 418.175e-03),
    (1.244818, 420.82e-03),
    (1.243191, 429.775e-03),
    (1.243773, 441.367e-03),
    (1.246324, 448.472e-03),
    (1.253665, 453.342e-03),
    (1.261358, 454.77e-03),
    (1.269539, 455.453e-03),
    (1.279789, 455.628e-03),
    (1.290039, 455.453e-03),
    (1.29804, 454.795e-03),
    (1.305914, 453.342e-03),
    (1.312847, 448.992e-03),
    (1.315806, 441.367e-03),
    (1.316355, 430.434e-03),
]

y_coords = [y for y, _ in rail_profile_points]
z_coords = [z for _, z in rail_profile_points]

y_center = (min(y_coords) + max(y_coords)) / 2
z_center = (min(z_coords) + max(z_coords)) / 2

rail_profile_points = [(y - y_center, z - z_center) for y, z in rail_profile_points]

doc = App.newDocument("Sweep3D")
profile_points = [App.Vector(0, y, z) for y, z in rail_profile_points]
profile_points.append(profile_points[0])
profile_wire_poly = Part.makePolygon(profile_points)
profile_face = Part.Face(profile_wire_poly)
profile_obj = doc.addObject("Part::Feature", "Profile")
profile_obj.Shape = profile_face
doc.recompute()

xs = []
ys = []
zs = []
us = []

txt_file = "C:/Users/bettif/Documents/Github/estra/out/results_0/models/Track.output/Track-Trk_Track.txt"
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

    for line in lines[xline + 3 : yline]:
        values = [float(value) for value in line.split(",")[:2]]
        xs.append(values[1])

    for line in lines[yline + 3 : zline]:
        values = [float(value) for value in line.split(",")[:2]]
        ys.append(values[1])

    for line in lines[zline + 3 : superelevation_line]:
        values = [float(value) for value in line.split(",")[:2]]
        zs.append(values[1])

    for line in lines[superelevation_line + 3 :]:
        values = [float(value) for value in line.split(",")[:2]]
        us.append(values[1])

points = [App.Vector(xs[i], ys[i], zs[i]) for i in range(len(xs))]

spline = Part.BSplineCurve()
spline.approximate(points)
edge = spline.toShape()
path_wire = Part.Wire([profile_wire_poly])
spline_obj = doc.addObject("Part::Feature", "Spine")
spline_obj.Shape = edge
doc.recompute()

# sweep = doc.addObject("Part::Sweep", "Sweep")
# sweep.Sections = [profile_obj]
# sweep.Spine = spline_obj
# sweep.Solid = True
# sweep.Frenet = False
# doc.recompute()

ps = Part.BRepOffsetAPI.MakePipeShell(path_wire)
spine_support = path_wire.extrude(App.Vector(0, 0, 1))
ps.setSpineSupport(spine_support)
ps.add(profile_wire_poly, False, True)

if ps.isReady():
    ps.build()
    ps.makeSolid()
    if ps.getStatus() == 0:
        tw = ps.shape()
        brep = Part.show(tw, "BREPPipeshell")


doc.recompute()

Gui.activeDocument().activeView().viewIsometric()
Gui.SendMsgToActiveView("ViewFit")

output_file = "C:/Users/bettif/Documents/Github/estra/out/results_0/models/Track.output/Track.step"
Part.export([brep], output_file)

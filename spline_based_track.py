from abaqus import *
from abaqusConstants import *
from assembly import *
from part import *
from material import *
from section import *
from step import *
from interaction import *
from load import *
from mesh import *
from optimization import *
from job import *
from sketch import *
from visualization import *
from connectorBehavior import *

import mesh
import regionToolset

import json
import os
import numpy as np
from scipy.interpolate import CubicSpline, interp1d

EPS = 1e-4


class AbaqusConstants:
    """Constants for Abaqus model element names."""

    MODEL_NAME = "SplineBasedTrack"
    RAIL_PART_NAME = "Rail"
    SLEEPER_PART_NAME = "Sleeper"

    CONCRETE_MATERIAL = "CONCRETE"
    STEEL_MATERIAL = "bbr1_ENG-Stainless Steel (ferritic)"

    RAIL_SECTION = "RailSection"
    SLEEPER_SECTION = "SleeperSection"

    TOP_RAIL_EDGE = "TopRailEdge"
    TOP_RAIL_NODES = "TopRailNodes"
    RAIL_TO_SLEEPER_COUPLING = "RailToSleeperCoupling"
    RAIL_TO_SLEEPER_COUPLING_CONTROL_POINT = "RailToSleeperCouplingControlPoint"
    LOWER_RAIL_EDGE = "LowerRailEdge"

    LEFT_BALLAST_CONNECTION = "LeftBallastConnection"
    RIGHT_BALLAST_CONNECTION = "RightBallastConnection"
    LEFT_BALLAST_CONNECTION_CONTROL_POINT = "LeftBallastConnectionControlPoint"
    RIGHT_BALLAST_CONNECTION_CONTROL_POINT = "RightBallastConnectionControlPoint"

    SLEEPER_TO_RAIL_COUPLING_LEFT = "SleeperToRailCouplingLeft"
    SLEEPER_TO_RAIL_COUPLING_RIGHT = "SleeperToRailCouplingRight"
    SLEEPER_TO_RAIL_COUPLING_LEFT_CONTROL_POINT = (
        "SleeperToRailCouplingLeftControlPoint"
    )
    SLEEPER_TO_RAIL_COUPLING_RIGHT_CONTROL_POINT = (
        "SleeperToRailCouplingRightControlPoint"
    )

    SLEEPER_TO_CONTROL_POINT_LEFT_COUPLING = "SleeperToControlPointLeftCoupling"
    SLEEPER_TO_CONTROL_POINT_RIGHT_COUPLING = "SleeperToControlPointRightCoupling"
    RAIL_TO_SLEEPER_COUPLING_PREFIX = "RailToSleeperCoupling"
    LEFT_BALLAST_CONNECTION_COUPLING = "LeftBallastConnectionCoupling"
    RIGHT_BALLAST_CONNECTION_COUPLING = "RightBallastConnectionCoupling"

    RAIL_TO_SLEEPER_LEFT_PREFIX = "RailToSleeperLeft"
    RAIL_TO_SLEEPER_RIGHT_PREFIX = "RailToSleeperRight"
    BALLAST_LEFT_PREFIX = "BallastLeft"
    BALLAST_RIGHT_PREFIX = "BallastRight"

    BEGINNING_RAIL_LOCKED = "BeginningRailLocked"
    ENDING_RAIL_LOCKED = "EndingRailLocked"

    BEGINNING_LEFT_RAIL_LOCKED = "BeginningLeftRailLocked"
    BEGINNING_RIGHT_RAIL_LOCKED = "BeginningRightRailLocked"
    ENDING_RIGHT_RAIL_LOCKED = "EndingRightRailLocked"
    ENDING_LEFT_RAIL_LOCKED = "EndingLeftRailLocked"

    INITIAL_STEP = "Initial"
    FREQUENCY_STEP = "FrequencyStep"

    SLEEPER_INSTANCE_PATTERN = "Sleeper-{}"
    RAIL_LEFT_INSTANCE = "Rail-Left"
    RAIL_RIGHT_INSTANCE = "Rail-Right"

    DATUM_CSYS_PATTERN = "DatumCsys{}"


class SplineBasedTrackGenerator:
    """Class to generate a spline-based track model from extracted geometry and spline data."""

    def __init__(self, model_name: str, json_file: str):
        """Initialize the spline-based track generator.

        Args:
            model_name (str): Name for the Abaqus model
        """
        self.model_name = model_name
        self.base_model_name = model_name
        self.model = None
        self.json_file = json_file
        with open(os.path.abspath(self.json_file), "r") as f:
            self.config = json.load(f)

        self.sleeper_length = self.config["SleeperLength"]
        self.rail_burn_in = self.config["RailBurnIn"]
        self.sleeper_spacing = self.config["SleeperSpacing"]
        self.rail_step_file = self.config["rail_step_file"]
        self.sleeper_step_file = self.config["sleeper_step_file"]
        self.cae_file_name = self.config["cae_file_name"]

        left_half_points = [
            (1.263218, 283.628e-03),
            (1.244588, 283.628e-03),
            (1.224584, 283.628e-03),
            (1.204789, 285.628e-03),
            (1.204789, 291.404e-03),
            (1.208504, 295.394e-03),
            (1.224584, 296.545e-03),
            (1.244588, 302.192e-03),
            (1.263218, 308.967e-03),
            (1.267464, 313.325e-03),
            (1.269523, 322.022e-03),
            (1.271539, 343.928e-03),
            (1.271539, 359.538e-03),
            (1.271131, 385.814e-03),
            (1.268339, 403.458e-03),
            (1.265682, 410.873e-03),
            (1.26189, 414.324e-03),
            (1.25177, 418.175e-03),
            (1.244818, 420.82e-03),
            (1.243191, 429.775e-03),
            (1.243773, 441.367e-03),
            (1.246324, 448.472e-03),
            (1.253665, 453.342e-03),
            (1.269539, 455.453e-03),
        ]

        center_top = [(1.279789, 455.628e-03)]
        center_bottom = [(1.279789, 283.628e-03)]

        def reflect_point(point):
            x, y = point
            center_x = center_top[0][0]
            new_x = 2 * center_x - x
            return (new_x, y)

        right_half_points = [reflect_point(p) for p in reversed(left_half_points)]
        rail_profile_points = (
            left_half_points + center_top + right_half_points + center_bottom
        )

        # Center the rail profile
        y_coords = [y for y, _ in rail_profile_points]
        z_coords = [z for _, z in rail_profile_points]
        y_center = (min(y_coords) + max(y_coords)) / 2
        z_center = (min(z_coords) + max(z_coords)) / 2
        self.rail_profile_points = [
            (y - y_center, z - z_center) for y, z in rail_profile_points
        ]
        self.rail_profile_points = [(-y, -z) for y, z in self.rail_profile_points]

        self.max_frequency = 200.0

        self.num_sleepers = 0
        self.lower_nodes = []
        self.txt_file = "Wheelset.output/Wheelset-Trk_Track.txt"
        self.segments = self.config.get(
            "segments",
            [
                {
                    "idx": 1,
                    "model_name": self.model_name,
                    "s_start": None,
                    "s_end": None,
                }
            ],
        )

    def get_spline_coordinates(
        self, file_path: str, s_start: float = None, s_end: float = None
    ):
        """Get spline coordinates from the track output file.

        Args:
            file_path (str): Path to the track output file

        Returns:
            tuple: Arrays of s, x, y, z coordinates
        """
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
                elif line.startswith('"Superelevation u(s)"'):
                    superelevation_line = idx

            for line in lines[xline + 3 : yline]:
                ss.append(float(line.split(",")[0]))
                xs.append(float(line.split(",")[1]))
            for line in lines[yline + 3 : zline]:
                ys.append(float(line.split(",")[1]))
            for line in lines[zline + 3 : superelevation_line]:
                zs.append(float(line.split(",")[1]))

        ss = np.array(ss)
        xs = np.array(xs)
        ys = np.array(ys)
        zs = np.array(zs)

        mask = np.ones_like(ss, dtype=bool)
        if s_start is not None:
            mask &= ss >= s_start
        if s_end is not None:
            mask &= ss <= s_end

        if not np.any(mask):
            raise ValueError(
                f"No spline data found between s_start={s_start} and s_end={s_end}"
            )

        ss = ss[mask]
        xs = xs[mask]
        ys = ys[mask]
        zs = zs[mask]

        ss = ss - ss[0]

        return ss, xs, ys, zs

    def find_top_point_for_wheel_contact(self, rail_part):
        """Find the top point of the rail profile for wheel contact.

        Returns:
            tuple: Top point coordinates (x, y, z)
        """
        min_x_coord = min([vertex.pointOn[0][0] for vertex in rail_part.vertices])
        zs = [z for _, z in self.rail_profile_points]
        top_z_index = zs.index(min(zs))
        top_point = (
            min_x_coord,
            self.rail_profile_points[top_z_index][0],
            self.rail_profile_points[top_z_index][1],
        )
        return top_point

    def find_lower_points_for_coupling(self, rail_part):
        """Find lower points of the rail profile for coupling.

        Returns:
            list: List of lower points
        """
        zs = [z for _, z in self.rail_profile_points]
        lower_points = []
        min_x_coord = min([vertex.pointOn[0][0] for vertex in rail_part.vertices])
        for i, (y, z) in enumerate(self.rail_profile_points):
            if np.isclose(z, max(zs)):
                lower_points.append((min_x_coord, y, z))
        return lower_points

    def find_edge_along_spline_given_initial_point(self, rail_part, initial_point):
        """Find edge along spline given initial point.

        Args:
            rail_part: Rail part object
            initial_point: Initial point coordinates

        Returns:
            int: Edge index
        """
        min_x_coord = min([vertex.pointOn[0][0] for vertex in rail_part.vertices])
        edge_to_find = None
        edges = rail_part.edges
        vertices = rail_part.vertices
        for edge in edges:
            vertices_indices = edge.getVertices()
            if (
                np.allclose(vertices[vertices_indices[0]].pointOn[0], initial_point)
                and vertices[vertices_indices[1]].pointOn[0][0] > min_x_coord
            ):
                edge_to_find = edge.index
                break
            elif (
                np.allclose(vertices[vertices_indices[1]].pointOn[0], initial_point)
                and vertices[vertices_indices[0]].pointOn[0][0] > min_x_coord
            ):
                edge_to_find = edge.index
                break
            else:
                continue
        return edge_to_find

    def create_model_and_parts(self, segment_idx: int):
        """Create the model and load parts from geometry files."""
        self.model = mdb.Model(name=self.model_name, modelType=STANDARD_EXPLICIT)
        mdb.openStep(self.rail_step_file, scaleFromFile=OFF)
        self.model.PartFromGeometryFile(
            combine=False,
            dimensionality=THREE_D,
            geometryFile=mdb.acis,
            name=AbaqusConstants.RAIL_PART_NAME,
            type=DEFORMABLE_BODY,
            bodyNum=segment_idx,
        )

        mdb.openStep(self.sleeper_step_file, scaleFromFile=OFF)
        self.model.PartFromGeometryFile(
            combine=False,
            dimensionality=THREE_D,
            geometryFile=mdb.acis,
            name=AbaqusConstants.SLEEPER_PART_NAME,
            type=DEFORMABLE_BODY,
            bodyNum=1,
        )

    def mesh_rail_part(self, mesh_size: float = 0.09):
        """Mesh the rail part with C3D8 elements.

        Args:
            mesh_size (float): Mesh size for the rail
        """
        rail_part = self.model.parts[AbaqusConstants.RAIL_PART_NAME]
        rail_part.seedPart(size=mesh_size, deviationFactor=0.5, minSizeFactor=0.5)

        rail_part.setElementType(
            elemTypes=(
                ElemType(
                    elemCode=C3D8,
                    elemLibrary=EXPLICIT,
                    secondOrderAccuracy=OFF,
                    distortionControl=DEFAULT,
                ),
            ),
            regions=(rail_part.cells,),
        )

        rail_part.generateMesh()

    def mesh_sleeper_part(self, mesh_size: float = 0.07):
        """Mesh the sleeper part with C3D6 elements.

        Args:
            mesh_size (float): Mesh size for the sleeper
        """
        sleeper_part = self.model.parts[AbaqusConstants.SLEEPER_PART_NAME]
        sleeper_part.seedPart(size=mesh_size, deviationFactor=0.5, minSizeFactor=0.5)

        sleeper_part.setElementType(
            elemTypes=(
                ElemType(
                    elemCode=C3D6,
                    elemLibrary=EXPLICIT,
                    secondOrderAccuracy=OFF,
                    distortionControl=DEFAULT,
                ),
            ),
            regions=(sleeper_part.cells,),
        )

        sleeper_part.generateMesh()

    def create_rail_sets_and_surfaces(self):
        """Create sets and surfaces for the rail part."""
        rail_part = self.model.parts[AbaqusConstants.RAIL_PART_NAME]

        top_point = self.find_top_point_for_wheel_contact(rail_part)
        lower_points = self.find_lower_points_for_coupling(rail_part)

        edges = rail_part.edges

        top_edge = self.find_edge_along_spline_given_initial_point(rail_part, top_point)
        edge_set = edges[top_edge : top_edge + 1]
        rail_part.Set(name=AbaqusConstants.TOP_RAIL_EDGE, edges=edge_set)
        nodes = rail_part.sets[AbaqusConstants.TOP_RAIL_EDGE].nodes
        rail_part.Set(name=AbaqusConstants.TOP_RAIL_NODES, nodes=nodes)
        del rail_part.sets[AbaqusConstants.TOP_RAIL_EDGE]

        self.lower_nodes = []
        for idx, lower_point in enumerate(lower_points):
            lower_edge = self.find_edge_along_spline_given_initial_point(
                rail_part, lower_point
            )
            edge_set = edges[lower_edge : lower_edge + 1]
            rail_part.Set(
                name=f"{AbaqusConstants.LOWER_RAIL_EDGE}{idx}", edges=edge_set
            )
            nodes = rail_part.sets[f"{AbaqusConstants.LOWER_RAIL_EDGE}{idx}"].nodes
            self.lower_nodes.extend(nodes)
            del rail_part.sets[f"{AbaqusConstants.LOWER_RAIL_EDGE}{idx}"]

    def add_coupling_sets_for_sleeper(self, x: float, j: int) -> bool:
        """Create rail-to-sleeper coupling sets based on sleeper positions.

        Args:
            ss: Spline parameter array
            xs_spline_interp: X coordinate spline interpolator
        """
        rail_part = self.model.parts[AbaqusConstants.RAIL_PART_NAME]

        labels_nodes_to_keep = []
        nodes_coordinates = []
        xs_start = x + self.rail_burn_in
        xs_end = x + self.sleeper_length + self.rail_burn_in

        for node in self.lower_nodes:
            if node.coordinates[0] >= xs_start and node.coordinates[0] <= xs_end:
                labels_nodes_to_keep.append(node.label)
                nodes_coordinates.append(node.coordinates)

        if len(labels_nodes_to_keep) > 0:
            nodes_to_keep = rail_part.nodes.sequenceFromLabels(
                labels=labels_nodes_to_keep
            )
            rail_part.Set(
                name=f"{AbaqusConstants.RAIL_TO_SLEEPER_COUPLING}{j+1}",
                nodes=nodes_to_keep,
            )

            control_point_label = rail_part.nodes.getClosest(
                coordinates=np.mean(nodes_coordinates, axis=0)
            ).label
            control_point = rail_part.nodes.sequenceFromLabels(
                labels=[control_point_label]
            )
            rail_part.Set(
                name=f"{AbaqusConstants.RAIL_TO_SLEEPER_COUPLING}{j+1}{AbaqusConstants.RAIL_TO_SLEEPER_COUPLING_CONTROL_POINT[-12:]}",
                nodes=control_point,
            )

        return len(labels_nodes_to_keep) > 0

    def create_sleeper_sets_and_surfaces(self):
        """Create sets and surfaces for the sleeper part."""
        sleeper_part = self.model.parts[AbaqusConstants.SLEEPER_PART_NAME]

        min_z_sleeper = max([node.coordinates[2] for node in sleeper_part.nodes])
        all_faces = sleeper_part.elementFaces

        element_faces = {"Left": [], "Right": []}
        faces_sides = {"Left": [], "Right": []}
        nodes_coordinates = {"Left": [], "Right": []}

        for face in all_faces:
            face_nodes = face.getNodes()
            nodes_z_coordinates = [node.coordinates[2] for node in face_nodes]

            if np.allclose(nodes_z_coordinates, min_z_sleeper) and np.allclose(face.getNormal(), [0., 0., 1.]):
                y_coordinates = face_nodes[0].coordinates[1]
                side = "Left" if y_coordinates < 0.0 else "Right"

                for node in face_nodes:
                    nodes_coordinates[side].append(node.coordinates)

                element_faces[side].append(face.getElements()[0])
                faces_sides[side].append(str(face.face))

        for side in ["Left", "Right"]:
            nodes_coordinates[side] = np.array(nodes_coordinates[side])
            control_point_label = sleeper_part.nodes.getClosest(
                coordinates=np.mean(nodes_coordinates[side], axis=0)
            ).label
            control_point = sleeper_part.nodes.sequenceFromLabels(
                labels=[control_point_label]
            )
            sleeper_part.Set(
                name=f"{side}{AbaqusConstants.LEFT_BALLAST_CONNECTION_CONTROL_POINT[4:]}",
                nodes=control_point,
            )

            surface_kwargs = {}
            for idx in range(1, 7):
                indices = [
                    index
                    for index, value in enumerate(faces_sides[side])
                    if value == f"FACE{idx}"
                ]
                elements = [element_faces[side][k] for k in indices]
                elements = sleeper_part.elements.sequenceFromLabels(
                    labels=[element.label for element in elements]
                )
                surface_kwargs[f"face{idx}Elements"] = elements

            surface_name = (
                AbaqusConstants.LEFT_BALLAST_CONNECTION
                if side == "Left"
                else AbaqusConstants.RIGHT_BALLAST_CONNECTION
            )
            sleeper_part.Surface(name=surface_name, **surface_kwargs)

        all_nodes = sleeper_part.nodes
        left_nodes = []
        right_nodes = []
        left_nodes_coordinates = []
        right_nodes_coordinates = []
        max_z_sleeper = min([node.coordinates[2] for node in all_nodes])

        for node in all_nodes:
            if np.isclose(node.coordinates[2], max_z_sleeper):
                if node.coordinates[1] < 0.0:
                    left_nodes.append(node)
                    left_nodes_coordinates.append(node.coordinates)
                else:
                    right_nodes.append(node)
                    right_nodes_coordinates.append(node.coordinates)

        left_nodes_sequence = sleeper_part.nodes.sequenceFromLabels(
            labels=[node.label for node in left_nodes]
        )
        sleeper_part.Set(
            name=AbaqusConstants.SLEEPER_TO_RAIL_COUPLING_LEFT,
            nodes=left_nodes_sequence,
        )

        right_nodes_sequence = sleeper_part.nodes.sequenceFromLabels(
            labels=[node.label for node in right_nodes]
        )
        sleeper_part.Set(
            name=AbaqusConstants.SLEEPER_TO_RAIL_COUPLING_RIGHT,
            nodes=right_nodes_sequence,
        )

        left_control_point_label = sleeper_part.nodes.getClosest(
            coordinates=np.mean(left_nodes_coordinates, axis=0)
        ).label
        left_control_point = sleeper_part.nodes.sequenceFromLabels(
            labels=[left_control_point_label]
        )
        sleeper_part.Set(
            name=AbaqusConstants.SLEEPER_TO_RAIL_COUPLING_LEFT_CONTROL_POINT,
            nodes=(left_control_point,),
        )

        right_control_point_label = sleeper_part.nodes.getClosest(
            coordinates=np.mean(right_nodes_coordinates, axis=0)
        ).label
        right_control_point = sleeper_part.nodes.sequenceFromLabels(
            labels=[right_control_point_label]
        )
        sleeper_part.Set(
            name=AbaqusConstants.SLEEPER_TO_RAIL_COUPLING_RIGHT_CONTROL_POINT,
            nodes=(right_control_point,),
        )

    def create_materials(self):
        """Create material properties."""
        concrete_material = self.model.Material(name=AbaqusConstants.CONCRETE_MATERIAL)
        concrete_material.Density(table=((2500.0,),))
        concrete_material.Elastic(table=((4.3e10, 0.15),))

        rail_material = self.model.Material(name=AbaqusConstants.STEEL_MATERIAL)
        rail_material.Density(table=((7800.0,),))
        rail_material.Conductivity(table=((18.0,),))
        rail_material.Elastic(table=((2e11, 0.28),))
        rail_material.Expansion(table=((1.1e-05,),))
        rail_material.Plastic(
            hardening=JOHNSON_COOK,
            table=((3.38674e7, 9.62742e8, 0.31129, 0.0, 0.0, 0.0),),
        )

    def create_sections(self):
        """Create sections and assign them to parts."""
        rail_section = self.model.HomogeneousSolidSection(
            name=AbaqusConstants.RAIL_SECTION,
            material=AbaqusConstants.STEEL_MATERIAL,
            thickness=None,
        )

        sleeper_section = self.model.HomogeneousSolidSection(
            name=AbaqusConstants.SLEEPER_SECTION,
            material=AbaqusConstants.CONCRETE_MATERIAL,
            thickness=None,
        )

        rail_region = regionToolset.Region(
            cells=self.model.parts[AbaqusConstants.RAIL_PART_NAME].cells
        )
        self.model.parts[AbaqusConstants.RAIL_PART_NAME].SectionAssignment(
            region=rail_region, sectionName=AbaqusConstants.RAIL_SECTION
        )

        sleeper_region = regionToolset.Region(
            cells=self.model.parts[AbaqusConstants.SLEEPER_PART_NAME].cells
        )
        self.model.parts[AbaqusConstants.SLEEPER_PART_NAME].SectionAssignment(
            region=sleeper_region, sectionName=AbaqusConstants.SLEEPER_SECTION
        )

    def create_assembly(self, ss, xs_spline_interp, ys_spline_interp, zs_spline_interp):
        """Create the assembly with sleepers and rails.

        Args:
            ss: Spline parameter array
            xs_spline_interp: X coordinate spline interpolator
            ys_spline_interp: Y coordinate spline interpolator
            zs_spline_interp: Z coordinate spline interpolator
        """
        assembly = self.model.rootAssembly

        # create rail instances
        rail_left_instance = assembly.Instance(
            name=AbaqusConstants.RAIL_LEFT_INSTANCE,
            part=self.model.parts[AbaqusConstants.RAIL_PART_NAME],
            dependent=ON,
        )
        rail_left_instance.translate(vector=(-self.rail_burn_in, -0.74, -0.22))

        rail_right_instance = assembly.Instance(
            name=AbaqusConstants.RAIL_RIGHT_INSTANCE,
            part=self.model.parts[AbaqusConstants.RAIL_PART_NAME],
            dependent=ON,
        )
        rail_right_instance.translate(vector=(-self.rail_burn_in, 0.77, -0.22))

        j = 0
        while (
            j * self.sleeper_spacing + self.sleeper_length + self.rail_burn_in
            < ss[-1] - EPS
        ):
            sleeper_instance = assembly.Instance(
                name=AbaqusConstants.SLEEPER_INSTANCE_PATTERN.format(j + 1),
                part=self.model.parts[AbaqusConstants.SLEEPER_PART_NAME],
                dependent=ON,
            )

            s = j * self.sleeper_spacing + EPS
            x, y, z = xs_spline_interp(s), ys_spline_interp(s), zs_spline_interp(s)

            deriv_y = (ys_spline_interp(s + EPS) - ys_spline_interp(s - EPS)) / (
                2 * EPS
            )
            deriv_x = (xs_spline_interp(s + EPS) - xs_spline_interp(s - EPS)) / (
                2 * EPS
            )

            tangent = np.array([deriv_x, deriv_y])
            tangent = tangent / np.linalg.norm(tangent)

            x, y, z = float(x), float(y), float(z)
            within_bounds = self.add_coupling_sets_for_sleeper(x, j)
            if within_bounds:
                sleeper_instance.translate(vector=(x, y, z))

                angle_about_y = np.rad2deg(np.arctan2(tangent[1], tangent[0]))
                sleeper_instance.rotateAboutAxis(
                    angle=angle_about_y,
                    axisDirection=(0, 0, 1),
                    axisPoint=(x, y, z),
                )

                self.num_sleepers += 1
                j += 1
            else:
                print(f"Broke after having added {self.num_sleepers} sleepers")
                break

    def create_coupling_constraints(self):
        """Create coupling constraints between rails and sleepers."""
        assembly = self.model.rootAssembly

        for i in range(self.num_sleepers):
            sleeper_instance = assembly.instances[
                AbaqusConstants.SLEEPER_INSTANCE_PATTERN.format(i + 1)
            ]
            left_rail_instance = assembly.instances[AbaqusConstants.RAIL_LEFT_INSTANCE]
            right_rail_instance = assembly.instances[
                AbaqusConstants.RAIL_RIGHT_INSTANCE
            ]

            # Sleeper to control point couplings
            self.model.Coupling(
                name=f"{sleeper_instance.name}.{AbaqusConstants.SLEEPER_TO_CONTROL_POINT_LEFT_COUPLING}",
                controlPoint=sleeper_instance.sets[
                    AbaqusConstants.SLEEPER_TO_RAIL_COUPLING_LEFT_CONTROL_POINT
                ],
                surface=sleeper_instance.sets[
                    AbaqusConstants.SLEEPER_TO_RAIL_COUPLING_LEFT
                ],
                influenceRadius=WHOLE_SURFACE,
                couplingType=DISTRIBUTING,
            )

            self.model.Coupling(
                name=f"{sleeper_instance.name}.{AbaqusConstants.SLEEPER_TO_CONTROL_POINT_RIGHT_COUPLING}",
                controlPoint=sleeper_instance.sets[
                    AbaqusConstants.SLEEPER_TO_RAIL_COUPLING_RIGHT_CONTROL_POINT
                ],
                surface=sleeper_instance.sets[
                    AbaqusConstants.SLEEPER_TO_RAIL_COUPLING_RIGHT
                ],
                influenceRadius=WHOLE_SURFACE,
                couplingType=DISTRIBUTING,
            )

            # Rail to sleeper couplings
            self.model.Coupling(
                name=f"{left_rail_instance.name}.LeftRailToSleeper{i+1}Coupling",
                controlPoint=left_rail_instance.sets[
                    f"{AbaqusConstants.RAIL_TO_SLEEPER_COUPLING}{i+1}{AbaqusConstants.RAIL_TO_SLEEPER_COUPLING_CONTROL_POINT[-12:]}"
                ],
                surface=left_rail_instance.sets[
                    f"{AbaqusConstants.RAIL_TO_SLEEPER_COUPLING}{i+1}"
                ],
                influenceRadius=WHOLE_SURFACE,
                couplingType=DISTRIBUTING,
            )

            self.model.Coupling(
                name=f"{right_rail_instance.name}.RightRailToSleeper{i+1}Coupling",
                controlPoint=right_rail_instance.sets[
                    f"{AbaqusConstants.RAIL_TO_SLEEPER_COUPLING}{i+1}{AbaqusConstants.RAIL_TO_SLEEPER_COUPLING_CONTROL_POINT[-12:]}"
                ],
                surface=right_rail_instance.sets[
                    f"{AbaqusConstants.RAIL_TO_SLEEPER_COUPLING}{i+1}"
                ],
                influenceRadius=WHOLE_SURFACE,
                couplingType=DISTRIBUTING,
            )

            # Ballast connection couplings
            self.model.Coupling(
                name=f"{sleeper_instance.name}.{AbaqusConstants.LEFT_BALLAST_CONNECTION_COUPLING}",
                controlPoint=sleeper_instance.sets[
                    AbaqusConstants.LEFT_BALLAST_CONNECTION_CONTROL_POINT
                ],
                surface=sleeper_instance.surfaces[
                    AbaqusConstants.LEFT_BALLAST_CONNECTION
                ],
                influenceRadius=WHOLE_SURFACE,
                couplingType=KINEMATIC,
                u1=ON,
                u2=ON,
                u3=ON,
                ur1=ON,
                ur2=ON,
                ur3=ON,
            )

            self.model.Coupling(
                name=f"{sleeper_instance.name}.{AbaqusConstants.RIGHT_BALLAST_CONNECTION_COUPLING}",
                controlPoint=sleeper_instance.sets[
                    AbaqusConstants.RIGHT_BALLAST_CONNECTION_CONTROL_POINT
                ],
                surface=sleeper_instance.surfaces[
                    AbaqusConstants.RIGHT_BALLAST_CONNECTION
                ],
                influenceRadius=WHOLE_SURFACE,
                couplingType=KINEMATIC,
                u1=ON,
                u2=ON,
                u3=ON,
                ur1=ON,
                ur2=ON,
                ur3=ON,
            )

    def create_connector_sections(self):
        """Create connector sections for rail-to-sleeper and sleeper-to-ballast connections."""
        rail_to_sleeper_values = [
            100000000000.0,
            100000000000.0,
            150000000.0,
            500000000000.0,
            500000000000.0,
            500000000000.0,
        ]
        sleeper_to_ballast_values = [
            500000000000.0,
            500000000000.0,
            300000000.0,
            500000000000.0,
            500000000000.0,
            500000000000.0,
        ]

        for i in range(self.num_sleepers):
            connector_sections = [
                f"{AbaqusConstants.RAIL_TO_SLEEPER_LEFT_PREFIX}{i+1}",
                f"{AbaqusConstants.RAIL_TO_SLEEPER_RIGHT_PREFIX}{i+1}",
                f"{AbaqusConstants.BALLAST_LEFT_PREFIX}{i+1}",
                f"{AbaqusConstants.BALLAST_RIGHT_PREFIX}{i+1}",
            ]

            for section_name in connector_sections:
                if (
                    "RailToSleeperLeft" in section_name
                    or "RailToSleeperRight" in section_name
                ):
                    elasticity_values = rail_to_sleeper_values
                else:
                    elasticity_values = sleeper_to_ballast_values

                behavior_options = []
                for j in range(6):
                    behavior_option = ConnectorElasticity(
                        behavior=LINEAR,
                        components=[j + 1],
                        coupling=UNCOUPLED,
                        dependencies=0,
                        frequencyDependency=OFF,
                        table=((elasticity_values[j],),),
                        temperatureDependency=OFF,
                    )
                    behavior_options.append(behavior_option)

                self.model.ConnectorSection(
                    name=section_name,
                    assembledType=BUSHING,
                    behaviorOptions=behavior_options,
                    defaultTolerance=ON,
                    extrapolation=CONSTANT,
                    integration=UNSPECIFIED,
                    materialFlowFactor=1.0,
                    regularization=0.03,
                    regularize=ON,
                    rotationalType=NONE,
                    translationalType=NONE,
                )

    def create_wire_polylines(self):
        """Create wire polylines for connector assignments."""
        assembly = self.model.rootAssembly

        for i in range(self.num_sleepers):
            sleeper_instance = assembly.instances[
                AbaqusConstants.SLEEPER_INSTANCE_PATTERN.format(i + 1)
            ]
            left_rail_instance = assembly.instances[AbaqusConstants.RAIL_LEFT_INSTANCE]
            right_rail_instance = assembly.instances[
                AbaqusConstants.RAIL_RIGHT_INSTANCE
            ]

            left_rail_control_point = left_rail_instance.sets[
                f"{AbaqusConstants.RAIL_TO_SLEEPER_COUPLING}{i+1}{AbaqusConstants.RAIL_TO_SLEEPER_COUPLING_CONTROL_POINT[-12:]}"
            ].nodes[0]
            left_sleeper_control_point = sleeper_instance.sets[
                AbaqusConstants.SLEEPER_TO_RAIL_COUPLING_LEFT_CONTROL_POINT
            ].nodes[0]
            assembly.WirePolyLine(
                points=((left_rail_control_point, left_sleeper_control_point),),
                mergeType=IMPRINT,
                meshable=False,
            )

            right_rail_control_point = right_rail_instance.sets[
                f"{AbaqusConstants.RAIL_TO_SLEEPER_COUPLING}{i+1}{AbaqusConstants.RAIL_TO_SLEEPER_COUPLING_CONTROL_POINT[-12:]}"
            ].nodes[0]
            right_sleeper_control_point = sleeper_instance.sets[
                AbaqusConstants.SLEEPER_TO_RAIL_COUPLING_RIGHT_CONTROL_POINT
            ].nodes[0]
            assembly.WirePolyLine(
                points=((right_rail_control_point, right_sleeper_control_point),),
                mergeType=IMPRINT,
                meshable=False,
            )

            left_ballast_control_point = sleeper_instance.sets[
                AbaqusConstants.LEFT_BALLAST_CONNECTION_CONTROL_POINT
            ].nodes[0]
            assembly.WirePolyLine(
                points=((left_ballast_control_point, None),),
                mergeType=IMPRINT,
                meshable=False,
            )

            right_ballast_control_point = sleeper_instance.sets[
                AbaqusConstants.RIGHT_BALLAST_CONNECTION_CONTROL_POINT
            ].nodes[0]
            assembly.WirePolyLine(
                points=((right_ballast_control_point, None),),
                mergeType=IMPRINT,
                meshable=False,
            )

            connector_sections = [
                f"{AbaqusConstants.RAIL_TO_SLEEPER_LEFT_PREFIX}{i+1}",
                f"{AbaqusConstants.RAIL_TO_SLEEPER_RIGHT_PREFIX}{i+1}",
                f"{AbaqusConstants.BALLAST_LEFT_PREFIX}{i+1}",
                f"{AbaqusConstants.BALLAST_RIGHT_PREFIX}{i+1}",
            ]

            for k in range(4):
                self.model.rootAssembly.DatumCsysByThreePoints(
                    name=AbaqusConstants.DATUM_CSYS_PATTERN.format(4 * i + k + 1),
                    coordSysType=CARTESIAN,
                    line1=(1.0, 0.0, 0.0),
                    line2=(0.0, 1.0, 0.0),
                    origin=(0.0, 0.0, 0.0),
                )

                index = self.model.rootAssembly.edges[k].index
                edge_array = self.model.rootAssembly.edges[index : index + 1]
                region = regionToolset.Region(edges=edge_array)

                datum_id = self.model.rootAssembly.datums.keys()[-1]
                section_assignment_name = connector_sections[k]

                self.model.rootAssembly.SectionAssignment(
                    region=region, sectionName=section_assignment_name
                )
                self.model.rootAssembly.ConnectorOrientation(
                    region=region,
                    axis1=AXIS_1,
                    axis2=AXIS_1,
                    angle1=0.0,
                    angle2=0.0,
                    localCsys1=self.model.rootAssembly.datums[datum_id],
                    orient2sameAs1=True,
                )

    def create_encastre_sets_for_rail(self):
        """Create encastre sets for the rail."""
        rail_part = self.model.parts[AbaqusConstants.RAIL_PART_NAME]
        nodes = rail_part.nodes

        node_labels = []
        for node in nodes:
            if np.isclose(node.coordinates[0], 0.0, atol=0.02):
                node_labels.append(node.label)

        nodes_to_be_locked = rail_part.nodes.sequenceFromLabels(labels=node_labels)
        rail_part.Set(
            name=AbaqusConstants.BEGINNING_RAIL_LOCKED, nodes=nodes_to_be_locked
        )

        max_x_coordinate = max([node.coordinates[0] for node in nodes])

        node_labels = []
        for node in nodes:
            if np.isclose(node.coordinates[0], max_x_coordinate, atol=0.02):
                node_labels.append(node.label)

        nodes_to_be_locked = rail_part.nodes.sequenceFromLabels(labels=node_labels)
        rail_part.Set(name=AbaqusConstants.ENDING_RAIL_LOCKED, nodes=nodes_to_be_locked)

    def create_boundary_conditions(self):
        """Create boundary conditions for the model."""
        assembly = self.model.rootAssembly

        left_rail_instance = assembly.instances[AbaqusConstants.RAIL_LEFT_INSTANCE]

        beginning_left_nodes = left_rail_instance.sets[
            AbaqusConstants.BEGINNING_RAIL_LOCKED
        ].nodes
        region = regionToolset.Region(nodes=beginning_left_nodes)
        self.model.EncastreBC(
            name=AbaqusConstants.BEGINNING_LEFT_RAIL_LOCKED,
            createStepName=AbaqusConstants.INITIAL_STEP,
            region=region,
        )

        end_left_nodes = left_rail_instance.sets[
            AbaqusConstants.ENDING_RAIL_LOCKED
        ].nodes
        region = regionToolset.Region(nodes=end_left_nodes)
        self.model.EncastreBC(
            name=AbaqusConstants.ENDING_LEFT_RAIL_LOCKED,
            createStepName=AbaqusConstants.INITIAL_STEP,
            region=region,
        )

        right_rail_instance = assembly.instances[AbaqusConstants.RAIL_RIGHT_INSTANCE]
        beginning_right_nodes = right_rail_instance.sets[
            AbaqusConstants.BEGINNING_RAIL_LOCKED
        ].nodes
        region = regionToolset.Region(nodes=beginning_right_nodes)
        self.model.EncastreBC(
            name=AbaqusConstants.BEGINNING_RIGHT_RAIL_LOCKED,
            createStepName=AbaqusConstants.INITIAL_STEP,
            region=region,
        )

        ending_right_nodes = right_rail_instance.sets[
            AbaqusConstants.ENDING_RAIL_LOCKED
        ].nodes
        region = regionToolset.Region(nodes=ending_right_nodes)
        self.model.EncastreBC(
            name=AbaqusConstants.ENDING_RIGHT_RAIL_LOCKED,
            createStepName=AbaqusConstants.INITIAL_STEP,
            region=region,
        )

    def create_frequency_step(self):
        """Create the frequency extraction step."""
        self.model.FrequencyStep(
            name=AbaqusConstants.FREQUENCY_STEP,
            previous=AbaqusConstants.INITIAL_STEP,
            maxEigen=self.max_frequency,
            eigensolver=AMS,
            normalization=MASS,
            acousticCoupling=AC_OFF,
        )

    def create_and_submit_job(self, job_name: str) -> None:
        """Create and submit the Abaqus job."""
        job = mdb.Job(
            name=job_name,
            model=self.model_name,
            type=ANALYSIS,
            resultsFormat=SIM,
        )

        job.writeInput()
        time.sleep(2)
        self.add_custom_line_to_inp(job_name)

        time.sleep(5)

        job = mdb.JobFromInputFile(job_name, f"{job_name}.inp")
        job.submit()
        job.waitForCompletion()

    def add_custom_line_to_inp(self, job_name: str) -> None:
        """Add a custom line to the generated .inp file."""
        inp_file_path = f"{job_name}.inp"
        
        with open(inp_file_path, 'r', encoding='utf-8') as f:
            lines = f.readlines()
            
        target_pattern = f", , {self.max_frequency:.0f}., , , ,"
        insert_position = None
        
        for i, line in enumerate(lines):
            if target_pattern in line:
                insert_position = i + 1
                break
        
        if insert_position is not None:
            custom_line = "*Flexible Body, type=SIMPACK\n"
            lines.insert(insert_position, custom_line)
            
            with open(inp_file_path, 'w', encoding='utf-8') as f:
                f.writelines(lines)

    def generate_model(self):
        """Generate the complete spline-based track models for each segment."""
        for segment in self.segments:
            segment_idx = segment.get("idx", 1)
            segment_name = segment.get(
                "model_name", f"{self.base_model_name}_Seg{segment_idx:02d}"
            )
            s_start = segment.get("s_start")
            s_end = segment.get("s_end")

            print(
                f"Generating segment {segment_idx}: {segment_name} "
                f"(s_start={s_start}, s_end={s_end})"
            )

            self.model_name = segment_name
            self.num_sleepers = 0
            self.lower_nodes = []

            ss, xs, ys, zs = self.get_spline_coordinates(
                self.txt_file, s_start, s_end
            )
            xs_spline_interp = interp1d(ss, xs)
            ys_spline_interp = interp1d(ss, ys)
            zs_spline_interp = interp1d(ss, zs)

            self.create_model_and_parts(segment_idx)

            self.mesh_rail_part()
            self.mesh_sleeper_part()

            self.create_rail_sets_and_surfaces()
            self.create_sleeper_sets_and_surfaces()

            self.create_materials()
            self.create_sections()

            self.create_assembly(ss, xs_spline_interp, ys_spline_interp, zs_spline_interp)

            self.create_coupling_constraints()
            self.create_connector_sections()
            self.create_wire_polylines()

            self.create_encastre_sets_for_rail()
            self.create_boundary_conditions()
            self.create_frequency_step()

            os.chdir(".")
            # self.create_and_submit_job(segment_name)

        mdb.saveAs(self.cae_file_name)


if __name__ == "__main__":
    generator = SplineBasedTrackGenerator(
        model_name=AbaqusConstants.MODEL_NAME,
        json_file="Wheelset.output\spline_track_config.json",
    )
    generator.generate_model()

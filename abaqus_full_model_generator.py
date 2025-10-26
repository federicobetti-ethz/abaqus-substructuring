"""Abaqus Script for Generating Single Sleeper Track Model"""

import numpy as np
import math
import random
import time

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
from scipy.interpolate import interp1d, CubicSpline

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


class AbaqusConstants:
    """Constants for Abaqus model element names."""
    
    MODEL_NAME = "SingleSleeperTrack"
    RAIL_PART_NAME = "Rail"
    SLEEPER_PART_NAME = "Sleeper"
    
    CONCRETE_MATERIAL = "CONCRETE"
    STEEL_MATERIAL = "bbr1_ENG-Stainless Steel (ferritic)"
    
    RAIL_SECTION = "RailSection"
    SLEEPER_SECTION = "SleeperSection"
    
    LEFT_BALLAST_CONNECTION = "LeftBallastConnection"
    RIGHT_BALLAST_CONNECTION = "RightBallastConnection"
    TIE_CONSTRAINT_BACKWARD = "TieConstraintBackward"
    TIE_CONSTRAINT_FORWARD = "TieConstraintForward"
    
    RAIL_TO_SLEEPER_COUPLING = "RailToSleeperCoupling"
    RAIL_TO_SLEEPER_COUPLING_CONTROL_POINT = "RailToSleeperCouplingControlPoint"
    SLEEPER_TO_RAIL_COUPLING_LEFT = "SleeperToRailCouplingLeft"
    SLEEPER_TO_RAIL_COUPLING_LEFT_CONTROL_POINT = "SleeperToRailCouplingLeftControlPoint"
    SLEEPER_TO_RAIL_COUPLING_RIGHT = "SleeperToRailCouplingRight"
    SLEEPER_TO_RAIL_COUPLING_RIGHT_CONTROL_POINT = "SleeperToRailCouplingRightControlPoint"
    LEFT_BALLAST_CONNECTION_CONTROL_POINT = "LeftBallastConnectionControlPoint"
    RIGHT_BALLAST_CONNECTION_CONTROL_POINT = "RightBallastConnectionControlPoint"
    
    SLEEPER_TO_CONTROL_POINT_LEFT_COUPLING = "SleeperToControlPointLeftCoupling"
    SLEEPER_TO_CONTROL_POINT_RIGHT_COUPLING = "SleeperToControlPointRightCoupling"
    RAIL_TO_CONTROL_POINT_COUPLING = "RailToControlPointCoupling"
    LEFT_BALLAST_CONNECTION_COUPLING = "LeftBallastConnectionCoupling"
    RIGHT_BALLAST_CONNECTION_COUPLING = "RightBallastConnectionCoupling"
    
    TIE_CONSTRAINT_PREFIX = "TieConstraint"
    
    RAIL_TO_SLEEPER_LEFT_PREFIX = "RailToSleeperLeft"
    RAIL_TO_SLEEPER_RIGHT_PREFIX = "RailToSleeperRight"
    BALLAST_LEFT_PREFIX = "BallastLeft"
    BALLAST_RIGHT_PREFIX = "BallastRight"
    RAIL_TO_WHEEL_PAIR_CONTACT = "RailToWheelPairContact"
    
    FIRST_LEFT_RAIL_LOCKED = "FirstLeftRailLocked"
    LAST_LEFT_RAIL_LOCKED = "LastLeftRailLocked"
    FIRST_RIGHT_RAIL_LOCKED = "FirstRightRailLocked"
    LAST_RIGHT_RAIL_LOCKED = "LastRightRailLocked"
    
    INITIAL_STEP = "Initial"
    FREQUENCY_STEP = "Frequency Step1"
    
    SLEEPER_INSTANCE_PATTERN = "Sleeper-{}"
    RAIL_LEFT_INSTANCE_PATTERN = "Rail-Left-{}"
    RAIL_RIGHT_INSTANCE_PATTERN = "Rail-Right-{}"
    
    DATUM_CSYS_PATTERN = "Datum-Csys-{}"
    
    CAE_FILE_NAME = "single_sleeper_track.cae"


class FlexTrackGenerator:
    """Class to generate a single sleeper track model from extracted data."""

    def __init__(
        self,
        model_name: str,
        input_json_data: dict,
        fe_data_folder: str,
    ):
        """Initialize the single sleeper generator.

        Args:
            model_name (str): Name for the Abaqus model
            input_json_data (dict): Input JSON data for the model
            fe_data_folder (str): Folder to save the FE data
        """
        self.model_name = model_name
        self.model = None
        self.input_json_data = input_json_data
        self.fe_data_folder = fe_data_folder

        self.num_repetitions = input_json_data["num_repetitions"]
        self.num_eigenmodes = input_json_data["num_eigenmodes"]

        self.sleeper_length = input_json_data["SleeperLength"]
        self.rail_burn_in = input_json_data["RailBurnIn"]
        self.rail_spacing = input_json_data["RailSpacing"]
        
        self.route_data = input_json_data["route"]

        self.sleeper_profile_points = sleeper_profile_points
        self.rail_profile_points = rail_profile_points

        self.rail_length = self.sleeper_length + 2 * self.rail_burn_in
        self.mesh_size_sleeper = input_json_data["meshSizeSleeper"]
        self.mesh_size_rail = input_json_data["meshSizeRail"]

        self.transform_matrix = [
            0.0, -1.0, 0.0,
            0.0, 0.0, -1.0,
            1.0, 0.0, 0.0,
            0.0, 0.0, 0.0,
        ]

    def create_model(self) -> mdb.Model:
        """Create a new Abaqus model."""
        if self.model_name in mdb.models.keys():
            del mdb.models[self.model_name]

        self.model = mdb.Model(name=self.model_name, modelType=STANDARD_EXPLICIT)
        return self.model

    def create_rail_profile_sketch(
        self,
        sketch_name: str
    ) -> object:
        """Create a sketch for the rail profile.

        Args:
            sketch_name (str): Name for the sketch

        Returns:
            object: Abaqus sketch object
        """
        transform_matrix = self.transform_matrix
        sketch = self.model.ConstrainedSketch(name=sketch_name, sheetSize=200.0, transform=transform_matrix)

        points = []
        for i, (y, z) in enumerate(self.rail_profile_points):
            points.append((y, z))

        for i in range(len(points) - 1):
            sketch.Line(point1=points[i], point2=points[i + 1])

        sketch.Line(point1=points[-1], point2=points[0])

        return sketch

    def create_sleeper_profile_sketch(
        self,
        sketch_name: str
    ) -> object:
        """
        Create a sketch for the ballast profile.

        Args:
            sketch_name (str): Name for the sketch

        Returns:
            object: Abaqus sketch object
        """
        transform_matrix = self.transform_matrix
        sketch = self.model.ConstrainedSketch(name=sketch_name, sheetSize=200.0, transform=transform_matrix)

        points = []
        for i, (y, z) in enumerate(self.sleeper_profile_points):
            points.append((y, z))

        for i in range(len(points) - 1):
            sketch.Line(point1=points[i], point2=points[i + 1])

        sketch.Line(point1=points[-1], point2=points[0])

        return sketch

    def create_rail_part(
        self,
        rail_name: str = AbaqusConstants.RAIL_PART_NAME,
    ) -> str:
        """Create a rail part by extruding the rail profile.

        Args:
            rail_name (str): Name for the rail part

        Returns:
            str: Name of the created rail part
        """
        rail_part = self.model.Part(name=rail_name, dimensionality=THREE_D, type=DEFORMABLE_BODY)
        sketch_name = f"rail_profile_{rail_name}"
        sketch = self.create_rail_profile_sketch(sketch_name)
        rail_part.BaseSolidExtrude(sketch=sketch, depth=self.rail_length)

        return rail_name

    def create_sleeper_part(
        self, sleeper_name: str = AbaqusConstants.SLEEPER_PART_NAME,
        sleeper_length: float = None,
    ) -> str:
        """Create a sleeper (ballast) part using the actual ballast profile.

        Args:
            sleeper_name (str): Name for the sleeper part

        Returns:
            str: Name of the created sleeper part
        """
        sleeper_part = self.model.Part(name=sleeper_name, dimensionality=THREE_D, type=DEFORMABLE_BODY)

        sketch_name = f"ballast_profile_{sleeper_name}"
        sketch = self.create_sleeper_profile_sketch(sketch_name)
        sleeper_part.BaseSolidExtrude(sketch=sketch, depth=self.sleeper_length)

        return sleeper_name

    def create_assembly(self,
        rail_part_name: str = AbaqusConstants.RAIL_PART_NAME,
        sleeper_part_name: str = AbaqusConstants.SLEEPER_PART_NAME,
        num_repetitions: int = None,
    ) -> None:
        """Create an assembly with the sleeper and two rails.

        Args:
            rail_part_name (str): Name of the rail part
            sleeper_part_name (str): Name of the sleeper part
            num_repetitions (int): Number of repetitions of the sleeper
        """
        assembly = self.model.rootAssembly

        for j in range(self.num_repetitions):
            sleeper_instance = assembly.Instance(
                name=AbaqusConstants.SLEEPER_INSTANCE_PATTERN.format(j + 1),
                part=self.model.parts[sleeper_part_name],
                dependent=ON,
            )
            rail1_instance = assembly.Instance(
                name=AbaqusConstants.RAIL_LEFT_INSTANCE_PATTERN.format(j + 1),
                part=self.model.parts[rail_part_name],
                dependent=ON,
            )
            rail2_instance = assembly.Instance(
                name=AbaqusConstants.RAIL_RIGHT_INSTANCE_PATTERN.format(j + 1),
                part=self.model.parts[rail_part_name],
                dependent=ON,
            )

            assembly.translate(
                instanceList=(AbaqusConstants.SLEEPER_INSTANCE_PATTERN.format(j + 1),),
                vector=(j * self.rail_length, 0, 0),
            )
            assembly.translate(
                instanceList=(AbaqusConstants.RAIL_LEFT_INSTANCE_PATTERN.format(j + 1),),
                vector=(-self.rail_burn_in + j * self.rail_length, 0, 0),
            )
            assembly.translate(
                instanceList=(AbaqusConstants.RAIL_RIGHT_INSTANCE_PATTERN.format(j + 1),),
                vector=(-self.rail_burn_in + j * self.rail_length, self.rail_spacing, 0),
            )

        for j in range(self.num_repetitions):
            self.create_coupling_constraints(
                AbaqusConstants.SLEEPER_INSTANCE_PATTERN.format(j + 1),
                AbaqusConstants.RAIL_LEFT_INSTANCE_PATTERN.format(j + 1),
                AbaqusConstants.RAIL_RIGHT_INSTANCE_PATTERN.format(j + 1),
            )

        for j in range(self.num_repetitions - 1):
            self.create_tie_constraints(
                AbaqusConstants.RAIL_LEFT_INSTANCE_PATTERN.format(j + 1),
                AbaqusConstants.RAIL_LEFT_INSTANCE_PATTERN.format(j + 2),
            )
            self.create_tie_constraints(
                AbaqusConstants.RAIL_RIGHT_INSTANCE_PATTERN.format(j + 1),
                AbaqusConstants.RAIL_RIGHT_INSTANCE_PATTERN.format(j + 2),
            )

        for j in range(self.num_repetitions):
            self.create_connector_sections_and_assignments(j + 1)
            self.create_wire_polylines_for_connectors(
                AbaqusConstants.SLEEPER_INSTANCE_PATTERN.format(j + 1),
                AbaqusConstants.RAIL_LEFT_INSTANCE_PATTERN.format(j + 1),
                AbaqusConstants.RAIL_RIGHT_INSTANCE_PATTERN.format(j + 1),
            )

        self.create_global_nodeset_for_rail_wheel_pair()
        self.create_boundary_conditions()
        # self.adjust_fe_model_to_route()

    def add_materials(self) -> None:
        """Add material properties for rails and sleeper using material parameters."""
        concrete_material = self.model.Material(name=AbaqusConstants.CONCRETE_MATERIAL)
        concrete_material.Density(table=((self.input_json_data["ConcreteDensity"],),))
        concrete_material.Elastic(
            table=(
                (
                    self.input_json_data["ConcreteYoungsModulus"],
                    self.input_json_data["ConcretePoissonRatio"],
                ),
            )
        )

        # Stainless Steel material (for rails)
        rail_material = self.model.Material(name=AbaqusConstants.STEEL_MATERIAL)
        rail_material.Density(table=((self.input_json_data["SteelDensity"],),))
        rail_material.Conductivity(table=((self.input_json_data["SteelConductivity"],),))
        rail_material.Elastic(
            table=(
                (
                    self.input_json_data["SteelYoungsModulus"],
                    self.input_json_data["SteelPoissonRatio"],
                ),
            )
        )
        rail_material.Expansion(table=((self.input_json_data["SteelExpansion"],),))
        rail_material.Plastic(
            hardening=JOHNSON_COOK,
            table=(
                (
                    self.input_json_data["SteelJCA"],
                    self.input_json_data["SteelJCB"],
                    self.input_json_data["SteelJCN"],
                    self.input_json_data["SteelJCM"],
                    self.input_json_data["SteelJCC"],
                    self.input_json_data["SteelJCEps0"],
                ),
            ),
        )

        # Create sections
        rail_section = self.model.HomogeneousSolidSection(
            name=AbaqusConstants.RAIL_SECTION, material=AbaqusConstants.STEEL_MATERIAL, thickness=None
        )

        sleeper_section = self.model.HomogeneousSolidSection(
            name=AbaqusConstants.SLEEPER_SECTION, material=AbaqusConstants.CONCRETE_MATERIAL, thickness=None
        )

        rail_region = regionToolset.Region(cells=self.model.parts[AbaqusConstants.RAIL_PART_NAME].cells)
        self.model.parts[AbaqusConstants.RAIL_PART_NAME].SectionAssignment(
            region=rail_region, sectionName=AbaqusConstants.RAIL_SECTION
        )

        sleeper_region = regionToolset.Region(cells=self.model.parts[AbaqusConstants.SLEEPER_PART_NAME].cells)
        self.model.parts[AbaqusConstants.SLEEPER_PART_NAME].SectionAssignment(
            region=sleeper_region, sectionName=AbaqusConstants.SLEEPER_SECTION
        )

    def mesh_rail_part(self, rail_part_name: str = AbaqusConstants.RAIL_PART_NAME) -> None:
        """Mesh the rail part with fine elements."""
        rail_part = self.model.parts[rail_part_name]
        rail_part.seedPart(size=self.mesh_size_rail, deviationFactor=0.1, minSizeFactor=0.2)
        rail_part.generateMesh()

        nodes = rail_part.nodes
        z_coordinates = [node.coordinates[2] for node in nodes]

        self.min_z_rail = min(z_coordinates)
        self.max_z_rail = max(z_coordinates)

    def mesh_sleeper_part(self, sleeper_part_name: str = AbaqusConstants.SLEEPER_PART_NAME) -> None:
        """Mesh the sleeper part with coarser elements."""
        sleeper_part = self.model.parts[sleeper_part_name]
        sleeper_part.seedPart(size=self.mesh_size_sleeper, deviationFactor=0.1, minSizeFactor=0.5)
        sleeper_part.generateMesh()

        nodes = sleeper_part.nodes
        z_coordinates = [node.coordinates[2] for node in nodes]

        self.min_z_sleeper = min(z_coordinates)
        self.max_z_sleeper = max(z_coordinates)

    def create_contact_nodeset_for_rail_wheel_pair(self, rail_part_name: str = AbaqusConstants.RAIL_PART_NAME) -> None:
        """Create a nodeset for the contact between the rail and wheel pair."""
        rail_part = self.model.parts[rail_part_name]
        all_nodes = rail_part.nodes

        nodes = []
        for node in all_nodes:
            if np.isclose(node.coordinates[2], self.min_z_rail):
                nodes.append(node)

        nodes = rail_part.nodes.sequenceFromLabels(labels=[node.label for node in nodes])
        rail_part.Set(name=AbaqusConstants.RAIL_TO_WHEEL_PAIR_CONTACT, nodes=nodes)

    def create_ballast_connection_surfaces(self, sleeper_part_name: str = AbaqusConstants.SLEEPER_PART_NAME) -> None:
        """Create surfaces for ballast connection on the bottom faces of the sleeper.

        Args:
            sleeper_part_name (str): Name of the sleeper part
        """
        sleeper_part = self.model.parts[sleeper_part_name]
        all_faces = sleeper_part.elementFaces

        element_faces = {"Left": [], "Right": []}
        faces_sides = {"Left": [], "Right": []}

        nodes_coordinates = {"Left": [], "Right": []}

        for face in all_faces:
            face_nodes = face.getNodes()
            nodes_z_coordinates = [node.coordinates[2] for node in face_nodes]

            if np.allclose(nodes_z_coordinates, self.max_z_sleeper):
                y_coordinates = face_nodes[0].coordinates[1]
                side = "Left" if y_coordinates < -0.5 else "Right"

                for node in face_nodes:
                    nodes_coordinates[side].append(node.coordinates)

                element_faces[side].append(face.getElements()[0])
                faces_sides[side].append(str(face.face))

        for side in ["Left", "Right"]:
            nodes_coordinates[side] = np.array(nodes_coordinates[side])
            control_point_label = sleeper_part.nodes.getClosest(
                coordinates=np.mean(nodes_coordinates[side], axis=0)
            ).label
            control_point = sleeper_part.nodes.sequenceFromLabels(labels=[control_point_label])
            sleeper_part.Set(name=f"{side}BallastConnectionControlPoint", nodes=control_point)

            face3Indices = [index for index, value in enumerate(faces_sides[side]) if value == "FACE3"]
            face3Elements = [element_faces[side][idx] for idx in face3Indices]
            face3Elements = sleeper_part.elements.sequenceFromLabels(
                labels=[element.label for element in face3Elements]
            )

            face4Indices = [index for index, value in enumerate(faces_sides[side]) if value == "FACE4"]
            face4Elements = [element_faces[side][idx] for idx in face4Indices]
            face4Elements = sleeper_part.elements.sequenceFromLabels(
                labels=[element.label for element in face4Elements]
            )

            face5Indices = [index for index, value in enumerate(faces_sides[side]) if value == "FACE5"]
            face5Elements = [element_faces[side][idx] for idx in face5Indices]
            face5Elements = sleeper_part.elements.sequenceFromLabels(
                labels=[element.label for element in face5Elements]
            )

            surface_kwargs = {}
            surface_kwargs["face3Elements"] = face3Elements
            surface_kwargs["face4Elements"] = face4Elements
            surface_kwargs["face5Elements"] = face5Elements

            surface_name = (
                AbaqusConstants.LEFT_BALLAST_CONNECTION if side == "Left" else AbaqusConstants.RIGHT_BALLAST_CONNECTION
            )
            sleeper_part.Surface(name=surface_name, **surface_kwargs)

    def create_tie_constraint_surfaces(self, rail_part_name: str = AbaqusConstants.RAIL_PART_NAME) -> None:
        """Create surfaces for tie constraints at the rail extrema (min and max y positions).

        Args:
            rail_part_name (str): Name of the rail part
        """
        rail_part = self.model.parts[rail_part_name]
        all_faces = rail_part.elementFaces

        backward_elements = []
        forward_elements = []

        backward_sides = []
        forward_sides = []

        for face in all_faces:
            face_nodes = face.getNodes()
            nodes_x_coordinates = [node.coordinates[0] for node in face_nodes]

            if np.allclose(nodes_x_coordinates, 0.0):
                backward_elements.append(face.getElements()[0])
                backward_sides.append(str(face.face))
            elif np.allclose(nodes_x_coordinates, self.rail_length):
                forward_elements.append(face.getElements()[0])
                forward_sides.append(str(face.face))

        backward_elements = rail_part.elements.sequenceFromLabels(
            labels=[element.label for element in backward_elements]
        )
        forward_elements = rail_part.elements.sequenceFromLabels(labels=[element.label for element in forward_elements])

        rail_part.Surface(name=AbaqusConstants.TIE_CONSTRAINT_BACKWARD, face1Elements=backward_elements)
        rail_part.Surface(name=AbaqusConstants.TIE_CONSTRAINT_FORWARD, face2Elements=forward_elements)

    def create_rail_coupling_surface(self, rail_part_name: str = AbaqusConstants.RAIL_PART_NAME) -> None:
        """Create a coupling surface and control point for the rail part.

        Args:
            rail_part_name (str): Name of the rail part
        """
        rail_part = self.model.parts[rail_part_name]
        all_nodes = rail_part.nodes

        nodes = []
        nodes_coordinates = []

        for node in all_nodes:
            if np.isclose(node.coordinates[2], self.max_z_rail):
                if np.logical_and(
                    node.coordinates[0] > self.rail_burn_in, node.coordinates[0] < self.rail_length - self.rail_burn_in
                ):
                    nodes.append(node)
                    nodes_coordinates.append(node.coordinates)

        nodes_coordinates = np.array(nodes_coordinates)

        nodes = rail_part.nodes.sequenceFromLabels(labels=[node.label for node in nodes])
        rail_part.Set(name=AbaqusConstants.RAIL_TO_SLEEPER_COUPLING, nodes=nodes)

        control_point_label = rail_part.nodes.getClosest(coordinates=np.mean(nodes_coordinates, axis=0)).label
        control_point = rail_part.nodes.sequenceFromLabels(labels=[control_point_label])
        rail_part.Set(name=AbaqusConstants.RAIL_TO_SLEEPER_COUPLING_CONTROL_POINT, nodes=control_point)

    def create_sleeper_coupling_surfaces(self, sleeper_part_name: str = AbaqusConstants.SLEEPER_PART_NAME) -> None:
        """Create coupling surfaces and control points for the left and right sides of the sleeper part.

        Args:
            sleeper_part_name (str): Name of the sleeper part
        """
        sleeper_part = self.model.parts[sleeper_part_name]
        all_nodes = sleeper_part.nodes

        left_nodes = []
        right_nodes = []

        left_nodes_coordinates = []
        right_nodes_coordinates = []

        for node in all_nodes:
            if np.isclose(node.coordinates[2], self.min_z_sleeper):
                if node.coordinates[1] < -0.5:
                    left_nodes.append(node)
                    left_nodes_coordinates.append(node.coordinates)
                else:
                    right_nodes.append(node)
                    right_nodes_coordinates.append(node.coordinates)

        left_nodes_sequence = sleeper_part.nodes.sequenceFromLabels(labels=[node.label for node in left_nodes])
        sleeper_part.Set(name=AbaqusConstants.SLEEPER_TO_RAIL_COUPLING_LEFT, nodes=left_nodes_sequence)

        right_nodes_sequence = sleeper_part.nodes.sequenceFromLabels(labels=[node.label for node in right_nodes])
        sleeper_part.Set(name=AbaqusConstants.SLEEPER_TO_RAIL_COUPLING_RIGHT, nodes=right_nodes_sequence)

        left_control_point_label = sleeper_part.nodes.getClosest(
            coordinates=np.mean(left_nodes_coordinates, axis=0)
        ).label
        left_control_point = sleeper_part.nodes.sequenceFromLabels(labels=[left_control_point_label])
        sleeper_part.Set(name=AbaqusConstants.SLEEPER_TO_RAIL_COUPLING_LEFT_CONTROL_POINT, nodes=(left_control_point,))

        right_control_point_label = sleeper_part.nodes.getClosest(
            coordinates=np.mean(right_nodes_coordinates, axis=0)
        ).label
        right_control_point = sleeper_part.nodes.sequenceFromLabels(labels=[right_control_point_label])
        sleeper_part.Set(
            name=AbaqusConstants.SLEEPER_TO_RAIL_COUPLING_RIGHT_CONTROL_POINT, nodes=(right_control_point,)
        )

    def create_coupling_constraints(
        self, sleeper_instance_name: str, left_rail_instance_name: str, right_rail_instance_name: str
    ) -> None:
        """Create coupling constraints between the rail and sleeper parts."""
        sleeper_instance = self.model.rootAssembly.instances[sleeper_instance_name]
        left_rail_instance = self.model.rootAssembly.instances[left_rail_instance_name]
        right_rail_instance = self.model.rootAssembly.instances[right_rail_instance_name]

        self.model.Coupling(
            name=f"{sleeper_instance_name}.{AbaqusConstants.SLEEPER_TO_CONTROL_POINT_LEFT_COUPLING}",
            controlPoint=sleeper_instance.sets[AbaqusConstants.SLEEPER_TO_RAIL_COUPLING_LEFT_CONTROL_POINT],
            surface=sleeper_instance.sets[AbaqusConstants.SLEEPER_TO_RAIL_COUPLING_LEFT],
            influenceRadius=WHOLE_SURFACE,
            couplingType=DISTRIBUTING,
        )

        self.model.Coupling(
            name=f"{sleeper_instance_name}.{AbaqusConstants.SLEEPER_TO_CONTROL_POINT_RIGHT_COUPLING}",
            controlPoint=sleeper_instance.sets[AbaqusConstants.SLEEPER_TO_RAIL_COUPLING_RIGHT_CONTROL_POINT],
            surface=sleeper_instance.sets[AbaqusConstants.SLEEPER_TO_RAIL_COUPLING_RIGHT],
            influenceRadius=WHOLE_SURFACE,
            couplingType=DISTRIBUTING,
        )

        self.model.Coupling(
            name=f"{left_rail_instance_name}.{AbaqusConstants.RAIL_TO_CONTROL_POINT_COUPLING}",
            controlPoint=left_rail_instance.sets[AbaqusConstants.RAIL_TO_SLEEPER_COUPLING_CONTROL_POINT],
            surface=left_rail_instance.sets[AbaqusConstants.RAIL_TO_SLEEPER_COUPLING],
            influenceRadius=WHOLE_SURFACE,
            couplingType=DISTRIBUTING,
        )

        self.model.Coupling(
            name=f"{right_rail_instance_name}.{AbaqusConstants.RAIL_TO_CONTROL_POINT_COUPLING}",
            controlPoint=right_rail_instance.sets[AbaqusConstants.RAIL_TO_SLEEPER_COUPLING_CONTROL_POINT],
            surface=right_rail_instance.sets[AbaqusConstants.RAIL_TO_SLEEPER_COUPLING],
            influenceRadius=WHOLE_SURFACE,
            couplingType=DISTRIBUTING,
        )

        self.model.Coupling(
            name=f"{sleeper_instance_name}.{AbaqusConstants.LEFT_BALLAST_CONNECTION_COUPLING}",
            controlPoint=sleeper_instance.sets[AbaqusConstants.LEFT_BALLAST_CONNECTION_CONTROL_POINT],
            surface=sleeper_instance.surfaces[AbaqusConstants.LEFT_BALLAST_CONNECTION],
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
            name=f"{sleeper_instance_name}.{AbaqusConstants.RIGHT_BALLAST_CONNECTION_COUPLING}",
            controlPoint=sleeper_instance.sets[AbaqusConstants.RIGHT_BALLAST_CONNECTION_CONTROL_POINT],
            surface=sleeper_instance.surfaces[AbaqusConstants.RIGHT_BALLAST_CONNECTION],
            influenceRadius=WHOLE_SURFACE,
            couplingType=KINEMATIC,
            u1=ON,
            u2=ON,
            u3=ON,
            ur1=ON,
            ur2=ON,
            ur3=ON,
        )

    def create_tie_constraints(self, backward_rail_instance_name: str, forward_rail_instance_name: str) -> None:
        """Create tie constraints between the backward and forward rail parts."""
        backward_rail_instance = self.model.rootAssembly.instances[backward_rail_instance_name]
        forward_rail_instance = self.model.rootAssembly.instances[forward_rail_instance_name]

        self.model.Tie(
            name=f"TieConstraint_{backward_rail_instance_name}_{forward_rail_instance_name}",
            main=backward_rail_instance.surfaces[AbaqusConstants.TIE_CONSTRAINT_FORWARD],
            secondary=forward_rail_instance.surfaces[AbaqusConstants.TIE_CONSTRAINT_BACKWARD],
            constraintEnforcement=SURFACE_TO_SURFACE,
            tieRotations=ON,
            positionToleranceMethod=COMPUTED,
            thickness=ON,
            adjust=OFF,
        )

    def create_wire_polylines_for_connectors(
        self, sleeper_instance_name: str, left_rail_instance_name: str, right_rail_instance_name: str
    ) -> None:
        """Create wire polylines for connector assignments."""
        assembly = self.model.rootAssembly

        sleeper_inst = assembly.instances[sleeper_instance_name]
        left_rail_inst = assembly.instances[left_rail_instance_name]
        right_rail_inst = assembly.instances[right_rail_instance_name]

        left_rail_control_point = left_rail_inst.sets[AbaqusConstants.RAIL_TO_SLEEPER_COUPLING_CONTROL_POINT].nodes[0]
        left_sleeper_control_point = sleeper_inst.sets[
            AbaqusConstants.SLEEPER_TO_RAIL_COUPLING_LEFT_CONTROL_POINT
        ].nodes[0]

        assembly.WirePolyLine(
            points=((left_rail_control_point, left_sleeper_control_point),), mergeType=IMPRINT, meshable=False
        )

        right_rail_control_point = right_rail_inst.sets[AbaqusConstants.RAIL_TO_SLEEPER_COUPLING_CONTROL_POINT].nodes[0]
        right_sleeper_control_point = sleeper_inst.sets[
            AbaqusConstants.SLEEPER_TO_RAIL_COUPLING_RIGHT_CONTROL_POINT
        ].nodes[0]

        assembly.WirePolyLine(
            points=((right_rail_control_point, right_sleeper_control_point),), mergeType=IMPRINT, meshable=False
        )

        left_ballast_control_point = sleeper_inst.sets[AbaqusConstants.LEFT_BALLAST_CONNECTION_CONTROL_POINT].nodes[0]

        assembly.WirePolyLine(points=((left_ballast_control_point, None),), mergeType=IMPRINT, meshable=False)

        right_ballast_control_point = sleeper_inst.sets[AbaqusConstants.RIGHT_BALLAST_CONNECTION_CONTROL_POINT].nodes[0]

        assembly.WirePolyLine(points=((right_ballast_control_point, None),), mergeType=IMPRINT, meshable=False)

        sleeper_number = int(sleeper_inst.name.split("-")[1])
        connector_sections = [
            f"{AbaqusConstants.RAIL_TO_SLEEPER_LEFT_PREFIX}_{sleeper_number}",
            f"{AbaqusConstants.RAIL_TO_SLEEPER_RIGHT_PREFIX}_{sleeper_number}",
            f"{AbaqusConstants.BALLAST_LEFT_PREFIX}_{sleeper_number}",
            f"{AbaqusConstants.BALLAST_RIGHT_PREFIX}_{sleeper_number}",
        ]

        for i in range(4):
            self.model.rootAssembly.DatumCsysByThreePoints(
                name=AbaqusConstants.DATUM_CSYS_PATTERN.format(4 * sleeper_number + i + 1),
                coordSysType=CARTESIAN,
                line1=(1.0, 0.0, 0.0),
                line2=(0.0, 1.0, 0.0),
                origin=(0.0, 0.0, 0.0),
            )

            index = self.model.rootAssembly.edges[i].index
            edge_array = self.model.rootAssembly.edges[index : index + 1]
            region = regionToolset.Region(edges=edge_array)

            datum_id = self.model.rootAssembly.datums.keys()[-1]
            section_assignment_name = connector_sections[i]

            self.model.rootAssembly.SectionAssignment(region=region, sectionName=section_assignment_name)

            self.model.rootAssembly.ConnectorOrientation(
                region=region,
                axis1=AXIS_1,
                axis2=AXIS_1,
                angle1=0.0,
                angle2=0.0,
                localCsys1=self.model.rootAssembly.datums[datum_id],
                orient2sameAs1=True,
            )

    def create_connector_sections_and_assignments(self, sleeper_number: int) -> None:
        """Create connector sections and assignments for the four specific wire polylines.

        Args:
            sleeper_number (int): Sleeper number for naming
        """
        connector_sections = [
            f"{AbaqusConstants.RAIL_TO_SLEEPER_LEFT_PREFIX}_{sleeper_number}",
            f"{AbaqusConstants.RAIL_TO_SLEEPER_RIGHT_PREFIX}_{sleeper_number}",
            f"{AbaqusConstants.BALLAST_LEFT_PREFIX}_{sleeper_number}",
            f"{AbaqusConstants.BALLAST_RIGHT_PREFIX}_{sleeper_number}",
        ]

        rail_to_sleeper_values = self.input_json_data["RailToSleeperElasticity"]
        sleeper_to_ballast_values = self.input_json_data["SleeperToBallastElasticity"]

        for section_name in connector_sections:
            if (
                AbaqusConstants.RAIL_TO_SLEEPER_LEFT_PREFIX in section_name
                or AbaqusConstants.RAIL_TO_SLEEPER_RIGHT_PREFIX in section_name
            ):
                elasticity_values = rail_to_sleeper_values
            else:
                elasticity_values = sleeper_to_ballast_values

            behavior_options = []
            for i in range(6):
                behavior_option = ConnectorElasticity(
                    behavior=LINEAR,
                    components=[i + 1],
                    coupling=UNCOUPLED,
                    dependencies=0,
                    frequencyDependency=OFF,
                    table=((elasticity_values[i],),),
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

    def create_global_nodeset_for_rail_wheel_pair(self) -> None:
        """Create a global nodeset for the contact between the rail and wheel pair."""
        right_sets = []
        for j in range(self.num_repetitions):
            right_rail_instance = self.model.rootAssembly.instances[AbaqusConstants.RAIL_RIGHT_INSTANCE_PATTERN.format(j + 1)]
            right_sets.append(right_rail_instance.sets[AbaqusConstants.RAIL_TO_WHEEL_PAIR_CONTACT])
        
        self.model.rootAssembly.SetByBoolean(name=f"Right{AbaqusConstants.RAIL_TO_WHEEL_PAIR_CONTACT}", sets=right_sets, operation=UNION)

        left_sets = []
        for j in range(self.num_repetitions):
            left_rail_instance = self.model.rootAssembly.instances[AbaqusConstants.RAIL_LEFT_INSTANCE_PATTERN.format(j + 1)]
            left_sets.append(left_rail_instance.sets[AbaqusConstants.RAIL_TO_WHEEL_PAIR_CONTACT])
        
        self.model.rootAssembly.SetByBoolean(name=f"Left{AbaqusConstants.RAIL_TO_WHEEL_PAIR_CONTACT}", sets=left_sets, operation=UNION)

    def create_boundary_conditions(self) -> None:
        """Create boundary conditions for the model."""
        first_left_rail_instance = self.model.rootAssembly.instances[
            AbaqusConstants.RAIL_LEFT_INSTANCE_PATTERN.format(1)
        ]
        last_left_rail_instance = self.model.rootAssembly.instances[
            AbaqusConstants.RAIL_LEFT_INSTANCE_PATTERN.format(self.num_repetitions)
        ]

        first_left_rail_locked_elements = first_left_rail_instance.surfaces[
            AbaqusConstants.TIE_CONSTRAINT_BACKWARD
        ].elements
        labels_nodes_to_be_locked = [
            node.label for element in first_left_rail_locked_elements for node in element.getNodes()
        ]
        labels_nodes_to_be_locked = list(dict.fromkeys(labels_nodes_to_be_locked))
        nodes_to_be_locked = first_left_rail_instance.nodes.sequenceFromLabels(labels=labels_nodes_to_be_locked)
        region = regionToolset.Region(nodes=nodes_to_be_locked)
        self.model.EncastreBC(
            name=AbaqusConstants.FIRST_LEFT_RAIL_LOCKED, createStepName=AbaqusConstants.INITIAL_STEP, region=region
        )

        last_left_rail_locked_elements = last_left_rail_instance.surfaces[
            AbaqusConstants.TIE_CONSTRAINT_FORWARD
        ].elements
        nodes_to_be_locked = [node.label for element in last_left_rail_locked_elements for node in element.getNodes()]
        nodes_to_be_locked = list(dict.fromkeys(nodes_to_be_locked))
        nodes_to_be_locked = last_left_rail_instance.nodes.sequenceFromLabels(labels=nodes_to_be_locked)
        region = regionToolset.Region(nodes=nodes_to_be_locked)
        self.model.EncastreBC(
            name=AbaqusConstants.LAST_LEFT_RAIL_LOCKED, createStepName=AbaqusConstants.INITIAL_STEP, region=region
        )

        first_right_rail_instance = self.model.rootAssembly.instances[
            AbaqusConstants.RAIL_RIGHT_INSTANCE_PATTERN.format(1)
        ]
        last_right_rail_instance = self.model.rootAssembly.instances[
            AbaqusConstants.RAIL_RIGHT_INSTANCE_PATTERN.format(self.num_repetitions)
        ]

        first_right_rail_locked_elements = first_right_rail_instance.surfaces[
            AbaqusConstants.TIE_CONSTRAINT_BACKWARD
        ].elements
        nodes_to_be_locked = [node.label for element in first_right_rail_locked_elements for node in element.getNodes()]
        nodes_to_be_locked = list(dict.fromkeys(nodes_to_be_locked))
        nodes_to_be_locked = first_right_rail_instance.nodes.sequenceFromLabels(labels=nodes_to_be_locked)
        region = regionToolset.Region(nodes=nodes_to_be_locked)
        self.model.EncastreBC(
            name=AbaqusConstants.FIRST_RIGHT_RAIL_LOCKED, createStepName=AbaqusConstants.INITIAL_STEP, region=region
        )

        last_right_rail_locked_elements = last_right_rail_instance.surfaces[
            AbaqusConstants.TIE_CONSTRAINT_FORWARD
        ].elements
        nodes_to_be_locked = [node.label for element in last_right_rail_locked_elements for node in element.getNodes()]
        nodes_to_be_locked = list(dict.fromkeys(nodes_to_be_locked))
        nodes_to_be_locked = last_right_rail_instance.nodes.sequenceFromLabels(labels=nodes_to_be_locked)
        region = regionToolset.Region(nodes=nodes_to_be_locked)
        self.model.EncastreBC(
            name=AbaqusConstants.LAST_RIGHT_RAIL_LOCKED, createStepName=AbaqusConstants.INITIAL_STEP, region=region
        )

    def adjust_fe_model_to_route(self) -> None:
        """Adjust the FE model to the route."""
        input_xs, input_elevations = self.route_data["x"], self.route_data["y"]
        interpolant = interp1d(input_xs, input_elevations, kind="linear", bounds_error=False, fill_value=(input_elevations[0], input_elevations[-1]))

        for j in range(self.num_repetitions):
            j_sleeper_position = j * self.rail_length
            slope = interpolant(j_sleeper_position) / j_sleeper_position if j_sleeper_position != 0 else 0.0
            rotation_angle = np.rad2deg(slope)

            sleeper_instance = self.model.rootAssembly.instances[AbaqusConstants.SLEEPER_INSTANCE_PATTERN.format(j + 1)]
            sleeper_instance.rotateAboutAxis(
                angle=rotation_angle,
                axisDirection=(0, 1, 0),
                axisPoint=(0, 0, 0),
            )

            left_rail_instance = self.model.rootAssembly.instances[AbaqusConstants.RAIL_LEFT_INSTANCE_PATTERN.format(j + 1)]
            left_rail_instance.rotateAboutAxis(
                angle=rotation_angle,
                axisDirection=(0, 1, 0),
                axisPoint=(0, 0, 0),
            )
            
            right_rail_instance = self.model.rootAssembly.instances[AbaqusConstants.RAIL_RIGHT_INSTANCE_PATTERN.format(j + 1)]
            right_rail_instance.rotateAboutAxis(
                angle=rotation_angle,
                axisDirection=(0, 1, 0),
                axisPoint=(0, 0, 0),
            )

    def define_frequency_extraction_step(self) -> None:
        """Define the frequency extraction step."""
        self.model.FrequencyStep(
            name=AbaqusConstants.FREQUENCY_STEP,
            previous=AbaqusConstants.INITIAL_STEP,
            numEigen=self.num_eigenmodes,
            eigensolver=LANCZOS,
            normalization=MASS,
            acousticCoupling=AC_PROJECTION,
        )

    def create_and_submit_job(self) -> None:
        """Create and submit the Abaqus job."""
        job = mdb.Job(
            name=self.model_name,
            model=self.model_name,
            type=ANALYSIS,
            resultsFormat=SIM,
        )

        job.writeInput()
        time.sleep(2)
        self.add_custom_line_to_inp()

        time.sleep(5)

        job = mdb.JobFromInputFile(self.model_name, f"{self.model_name}.inp")
        job.submit()
        job.waitForCompletion()

    def add_custom_line_to_inp(self) -> None:
        """Add a custom line to the generated .inp file."""
        inp_file_path = f"{self.model_name}.inp"
        
        with open(inp_file_path, 'r') as f:
            lines = f.readlines()
            
        target_pattern = f"{self.num_eigenmodes}, , , , , "
        insert_position = None
        
        for i, line in enumerate(lines):
            if target_pattern in line:
                insert_position = i + 1
                break
        
        if insert_position is not None:
            custom_line = "*Flexible Body, type=SIMPACK\n"
            lines.insert(insert_position, custom_line)
            
            with open(inp_file_path, 'w') as f:
                f.writelines(lines)

    def generate_model(self) -> None:
        """Generate the complete single sleeper track model."""
        self.create_model()
        self.create_rail_part()
        self.create_sleeper_part()
        
        self.add_materials()
        self.mesh_rail_part()
        self.mesh_sleeper_part()

        self.create_ballast_connection_surfaces()
        self.create_tie_constraint_surfaces()
        self.create_rail_coupling_surface()
        self.create_sleeper_coupling_surfaces()
        self.create_contact_nodeset_for_rail_wheel_pair()

        self.create_assembly()
        self.define_frequency_extraction_step()

        os.chdir(self.fe_data_folder)
        mdb.saveAs(pathName=self.model_name + ".cae")
        
        if self.input_json_data["solve"]:
            self.create_and_submit_job()

def main(model_name: str, input_json_data: dict, fe_data_folder: str):
    """Main function to generate the single sleeper track model."""
    generator = FlexTrackGenerator(
        model_name, 
        input_json_data=input_json_data, 
        fe_data_folder=fe_data_folder
    )
    generator.generate_model()

if __name__ == "__main__":
    model_name = "SingleSleeperTrack"
    fe_data_folder = r"out/results_0\fe-data"
    input_json_data = r"out/results_0\FlexTrack.json"

    with open(input_json_data, 'r') as f:
        input_json_data = json.load(f)
    
    main(
        model_name=model_name, 
        input_json_data=input_json_data, 
        fe_data_folder=fe_data_folder
    )

"""Example of how to use the FlexTrack class to create an Abaqus model."""

from estra.dynamics.abaqus.abaqus_track_model import AbaqusTrackModel
from estra.dynamics.abaqus.material_parameters import MaterialParameters
from estra.dynamics.abaqus.geometry_parameters import GeometryParameters

from estra import TEST_ROUTE_DATA
from estra.dynamics.simpack.utils.constants import MAX_DISTANCE_TRAVEL
from estra.datatypes.routes.route_from_csv import RouteFromCSV
from estra.dynamics.abaqus.abaqus_track_model import AbaqusTrackModel
from estra.dynamics.abaqus.material_parameters import MaterialParameters
from estra.dynamics.abaqus.geometry_parameters import GeometryParameters
import numpy as np
import torch

from pathlib import Path
from estra.datatypes.routes.track_spline_route import TrackSplineRoute

print("=" * 80)
print("ESTRA ABAQUS MODEL EXAMPLE")
print("=" * 80)

track_data_path = Path("out/results_0/models/Train.output/Train-Trk_Track.txt")

route = TrackSplineRoute(
    name="Train Track Spline",
    track_data_path=track_data_path,
    interpolation_method="cubic",
)

material_parameters = MaterialParameters()
geometry_parameters = GeometryParameters()
abaqus_track_model = AbaqusTrackModel(
    model_name="SingleSleeperTrack",
    material_parameters=material_parameters,
    geometry_parameters=geometry_parameters,
    route=route,
    length=100,
    num_eigenmodes=25,
    solve=False,
)
abaqus_track_model.create_abaqus_model()

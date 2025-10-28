from abaqus import *
from abaqusConstants import *
from interaction import *
from material import *
from part import *
from assembly import *
from step import *
import regionToolset as rts
from connectorBehavior import *

import time


def BASELINE_MAP_COUPLINGS(i):
    """Mapping of original model constraint names to new model constraint names."""
    range_to_consider = range(4 * (20 - i) + 1, 4 * (20 - i + 1) + 1)

    return {
        range_to_consider[0]: [
            "Coupling 1-1-{}-4".format(i),
            "Coupling 1-2-1-{}-1".format(i),
        ],
        range_to_consider[1]: [
            "Coupling 1-1-{}-1".format(i),
            "Coupling 1-1-{}-2".format(i),
        ],
        range_to_consider[2]: ["Coupling 1-1-{}".format(i)],
        range_to_consider[3]: ["Coupling 1-2-1-{}".format(i)],
    }


MAP_COUPLINGS = BASELINE_MAP_COUPLINGS(20)
for i in reversed(range(1, 20)):
    MAP_COUPLINGS.update(BASELINE_MAP_COUPLINGS(i))


def open_source_model(path, source_model_name):
    """Open the source model from the given path and return the model object.

    Args:
        path: Path to the source model.
        source_model_name: Name of the source model.

    Returns:
        Model object.
    """
    mdb = openMdb(path)
    src_model = mdb.models[source_model_name]
    return src_model


def add_tying_surfaces_for_replication(src_model, part_name, num_copies):
    """Add surfaces for boundary conditions to the original model.

    Args:
        src_model: Model object.
        part_name: Name of the part.
        num_copies: Number of copies.
    """
    # on the original model, add surfaces for boundary conditions (used for Tie couplings between track pieces)
    base_instance = src_model.rootAssembly.instances["{}-1".format(part_name)]
    if num_copies > 1:
        for bc_name, bc_obj in src_model.boundaryConditions.items():
            region_set = src_model.rootAssembly.sets[bc_obj.region[0]]
            region_elements_labels = [
                elem.label for node in region_set.nodes for elem in node.getElements()
            ]
            region_elements_labels = list(
                dict.fromkeys(region_elements_labels)
            )  # remove duplicates
            region_elements = base_instance.elements.sequenceFromLabels(
                region_elements_labels
            )
            region_elements_mask = region_elements.getMask()
            params = {"name": "{}".format(bc_name)}
            if str(bc_name)[-1] in ["1", "2"]:
                params["face1Elements"] = base_instance.elements.getSequenceFromMask(
                    mask=region_elements_mask
                )
            elif str(bc_name)[-1] in ["3", "4"]:
                params["face2Elements"] = base_instance.elements.getSequenceFromMask(
                    mask=region_elements_mask
                )
            src_model.rootAssembly.Surface(**params)


def define_new_model(new_model_name, part_name, src_part):
    """Define a new model with the given name and part.

    Args:
        new_model_name: Name of the new model.
        part_name: Name of the part.
        src_part: Part object.
    """
    new_model = mdb.Model(name=new_model_name)
    new_model.Part(name=part_name, objectToCopy=src_part)
    return new_model


def copy_part_and_translate(new_asm, part_name, part_copy, track_length, num_copies):
    """Copy the part and translate it.

    Args:
        new_asm: Assembly object.
        part_name: Name of the part.
        part_copy: Part object.
        track_length: Length of the track.
        num_copies: Number of copies.
    """
    for i in range(num_copies):
        inst_name = "{}-{}".format(part_name, i + 1)
        new_asm.Instance(name=inst_name, part=part_copy, dependent=ON)
        new_asm.translate(instanceList=(inst_name,), vector=(-i * track_length, 0, 0))


def define_model_materials(new_model):
    """Define the materials for the model. #TODO: take the values from the original model instead of hardcoding them.

    Args:
        new_model: Model object.
    """
    new_model.Material(name="CONCRETE")
    new_model.materials["CONCRETE"].Density(table=((2500.0,),))
    new_model.materials["CONCRETE"].Elastic(table=((4.3e10, 0.15),))

    new_model.Material(name="bbr1_ENG-Stainless Steel (ferritic)")
    new_model.materials["bbr1_ENG-Stainless Steel (ferritic)"].Density(
        table=((7800.0,),)
    )
    new_model.materials["bbr1_ENG-Stainless Steel (ferritic)"].Conductivity(
        table=((18.0,),)
    )
    new_model.materials["bbr1_ENG-Stainless Steel (ferritic)"].Elastic(
        table=((2e11, 0.28),)
    )
    new_model.materials["bbr1_ENG-Stainless Steel (ferritic)"].Expansion(
        table=((1.1e-5,),)
    )

    new_model.materials["bbr1_ENG-Stainless Steel (ferritic)"].Plastic(
        hardening=JOHNSON_COOK, table=((3.38674e7, 9.62742e8, 0.31129, 0.0, 0.0, 0.0),)
    )


def copy_sections_from_original_model(src_model, new_model, num_copies, part_name):
    """Copy the sections from the original model to the new model.

    Args:
        src_model: Model object.
        new_model: Model object.
        num_copies: Number of copies.
        part_name: Name of the part.
    """
    for section_name, section_obj in src_model.sections.items():
        if section_name.startswith("Sec"):
            new_model.HomogeneousSolidSection(
                name=section_name,
                material=section_obj.material,
                thickness=section_obj.thickness,
            )
        elif section_name.startswith("Con"):
            behaviorOptions = []
            for i in range(6):
                behavior_options_to_copy = section_obj.behaviorOptions[i]
                behavior_options = ConnectorElasticity(
                    behavior=behavior_options_to_copy.behavior,
                    components=[i + 1],
                    coupling=behavior_options_to_copy.coupling,
                    dependencies=behavior_options_to_copy.dependencies,
                    frequencyDependency=behavior_options_to_copy.frequencyDependency,
                    table=behavior_options_to_copy.table,
                    temperatureDependency=behavior_options_to_copy.temperatureDependency,
                )
                behaviorOptions.append(behavior_options)

            new_model.ConnectorSection(
                name=section_name,
                assembledType=BUSHING,
                behaviorOptions=behaviorOptions,
                defaultTolerance=ON,
                extrapolation=CONSTANT,
                integration=UNSPECIFIED,
                materialFlowFactor=1.0,
                regularization=0.03,
                regularize=ON,
                rotationalType=NONE,
                translationalType=NONE,
            )


def copy_sets_and_internal_sets_from_original_model(
    src_asm, new_asm, num_copies, part_name
):
    """Copy the sets and internal sets from the original model to the new model.

    Args:
        src_asm: Assembly object.
        new_asm: Assembly object.
        num_copies: Number of copies.
        part_name: Name of the part.
    """
    # copying all the sets to the new model
    for set_name, aset in src_asm.sets.items():
        set_mask = aset.nodes.getMask()
        for i in range(num_copies):
            inst_name = "{}-{}".format(part_name, i + 1)
            if set_mask:
                new_asm.Set(
                    name="{}_{}".format(inst_name, set_name),
                    nodes=new_asm.instances[inst_name].nodes.getSequenceFromMask(
                        mask=set_mask
                    ),
                )

    # all the NODE internal sets to the new model
    for internal_set_name, internal_aset in src_asm.allInternalSets.items():
        for i in range(num_copies):
            inst_name = "{}-{}".format(part_name, i + 1)
            if (
                internal_set_name.startswith("RAI")
                or internal_set_name.startswith("SL_")
                or internal_set_name.startswith("_M")
            ):
                node_set_mask = internal_aset.nodes.getMask()
                new_asm.Set(
                    name="{}_{}".format(inst_name, internal_set_name),
                    nodes=new_asm.instances[inst_name].nodes.getSequenceFromMask(
                        mask=node_set_mask
                    ),
                    internal=ON,
                )

    for internal_set_name, internal_aset in src_asm.allInternalSets.items():
        element_set_mask = internal_aset.elements.getMask()
        for i in range(num_copies):
            inst_name = "{}-{}".format(part_name, i + 1)
            if internal_set_name.startswith("_Picked") and element_set_mask:
                new_asm.Set(
                    name="{}_{}".format(inst_name, internal_set_name),
                    elements=new_asm.instances[inst_name].elements.getSequenceFromMask(
                        mask=element_set_mask
                    ),
                    internal=ON,
                )


def copy_surfaces_from_original_model(
    src_asm, new_asm, src_model, num_copies, part_name
):
    """Copy the surfaces from the original model to the new model.

    Args:
        src_asm: Assembly object.
        new_asm: Assembly object.
        src_model: Model object.
        num_copies: Number of copies.
        part_name: Name of the part.
    """
    # copying all the surfaces to the new model
    for surf_name, asurf in src_asm.surfaces.items():
        face_elements = asurf.elements
        if surf_name.startswith("TIE"):  # tie couplings between rail sections
            for i in range(num_copies):
                inst_name = "{}-{}".format(part_name, i + 1)
                elements_to_tag = new_asm.instances[
                    inst_name
                ].elements.getSequenceFromMask(mask=face_elements.getMask())
                if str(asurf.sides[0]) == "FACE1":
                    new_asm.Surface(
                        name="{}_{}".format(inst_name, surf_name),
                        face1Elements=elements_to_tag,
                    )
                elif str(asurf.sides[0]) == "FACE2":
                    new_asm.Surface(
                        name="{}_{}".format(inst_name, surf_name),
                        face2Elements=elements_to_tag,
                    )
                elif str(asurf.sides[0]) == "FACE4":
                    new_asm.Surface(
                        name="{}_{}".format(inst_name, surf_name),
                        face4Elements=elements_to_tag,
                    )

        elif surf_name.startswith("BAL"):  # ballast connections
            str_asurf_sides = [str(side) for side in asurf.sides]
            face4_tags = [
                index
                for index, element in enumerate(str_asurf_sides)
                if element == "FACE4"
            ]
            face5_tags = [
                index
                for index, element in enumerate(str_asurf_sides)
                if element == "FACE5"
            ]
            face4_mask = face_elements[face4_tags[0] : face5_tags[-1] + 1].getMask()
            face5_mask = face_elements[face5_tags[0] : face5_tags[-1] + 1].getMask()
            for i in range(num_copies):
                inst_name = "{}-{}".format(part_name, i + 1)
                elements_to_tag_face4 = new_asm.instances[
                    inst_name
                ].elements.getSequenceFromMask(mask=face4_mask)
                elements_to_tag_face5 = new_asm.instances[
                    inst_name
                ].elements.getSequenceFromMask(mask=face5_mask)
                new_asm.Surface(
                    name="{}_{}".format(inst_name, surf_name),
                    face4Elements=elements_to_tag_face4,
                    face5Elements=elements_to_tag_face5,
                )

        elif surf_name.startswith("Disp"):  # displacement boundary conditions
            for i in range(num_copies):
                inst_name = "{}-{}".format(part_name, i + 1)
                inst_instance = new_asm.instances[inst_name]
                surface_name = "{}_{}".format(inst_name, surf_name)
                base_params = {"name": surface_name}
                face_mask = src_model.rootAssembly.surfaces[
                    surf_name
                ].elements.getMask()
                if str(surf_name)[-1] in ["1", "2"]:
                    base_params["face1Elements"] = (
                        inst_instance.elements.getSequenceFromMask(mask=face_mask)
                    )
                elif str(surf_name)[-1] in ["3", "4"]:
                    base_params["face2Elements"] = (
                        inst_instance.elements.getSequenceFromMask(mask=face_mask)
                    )
                new_asm.Surface(**base_params)


def define_and_remove_redundant_bcs(src_model, new_model, num_copies, part_name):
    """Define and remove redundant boundary conditions.

    Args:
        src_model: Model object.
        new_model: Model object.
        num_copies: Number of copies.
        part_name: Name of the part.
    """
    # define all possible boundary conditions (more than needed, the ones that are not needed are deleted afterwards)
    for bc_name, bc_obj in src_model.boundaryConditions.items():
        bc_nodes_mask = src_model.rootAssembly.sets[bc_obj.region[0]].nodes.getMask()
        for i in range(num_copies):
            inst_name = "{}-{}".format(part_name, i + 1)
            inst_instance = new_model.rootAssembly.instances[inst_name]
            region = rts.Region(
                nodes=inst_instance.nodes.getSequenceFromMask(mask=bc_nodes_mask)
            )
            new_model.EncastreBC(
                name="{}_{}".format(inst_name, bc_name),
                createStepName="Initial",
                region=region,
            )

    # delete all boundary conditions except Disp-BC-3 and Disp-BC-4 for the first instance and Disp-BC-1 and Disp-BC-2 for the last instance
    if num_copies > 1:
        for i in range(num_copies):
            inst_name = "{}-{}".format(part_name, i + 1)
            if i == 0:
                bc_to_delete = ["Disp-BC-3", "Disp-BC-4"]
            elif i == num_copies - 1:
                bc_to_delete = ["Disp-BC-1", "Disp-BC-2"]
            else:
                bc_to_delete = ["Disp-BC-1", "Disp-BC-2", "Disp-BC-3", "Disp-BC-4"]
            for bc_name in bc_to_delete:
                del new_model.boundaryConditions["{}_{}".format(inst_name, bc_name)]


def copy_constraints_from_original_model(src_model, new_model, num_copies, part_name):
    """Copy the constraints from the original model to the new model.

    Args:
        src_model: Model object.
        new_model: Model object.
        num_copies: Number of copies.
        part_name: Name of the part.
    """
    for constr_name, constr_obj in src_model.constraints.items():
        if constr_name.startswith("Cou"):
            control_point_mask = src_model.rootAssembly.sets[
                constr_obj.controlPoint[0]
            ].nodes.getMask()
            surface_string = constr_obj.surface[0]
            for i in range(num_copies):
                inst_name = "{}-{}".format(part_name, i + 1)
                control_point_region = rts.Region(
                    nodes=new_model.rootAssembly.instances[
                        inst_name
                    ].nodes.getSequenceFromMask(mask=control_point_mask)
                )
                base_coupling_params = {
                    "name": "{}_{}".format(inst_name, constr_name),
                    "controlPoint": control_point_region,
                    "couplingType": constr_obj.couplingType,
                    "adjust": False,
                    "weightingMethod": UNIFORM,
                    "influenceRadius": WHOLE_SURFACE,
                }
                if str(constr_obj.couplingType) == "KINEMATIC":
                    base_coupling_params.update(
                        {
                            "u1": ON,
                            "u2": ON,
                            "u3": ON,
                            "ur1": ON,
                            "ur2": ON,
                            "ur3": ON,
                        }
                    )
                try:
                    base_coupling_params["surface"] = new_model.rootAssembly.surfaces[
                        "{}_{}".format(inst_name, surface_string)
                    ]
                    new_model.Coupling(**base_coupling_params)
                except:
                    base_coupling_params["surface"] = new_model.rootAssembly.sets[
                        "{}_{}".format(inst_name, surface_string)
                    ]
                    new_model.Coupling(**base_coupling_params)
        elif constr_name.startswith("TIE"):
            for i in range(num_copies):
                inst_name = "{}-{}".format(part_name, i + 1)
                new_model.Tie(
                    name="{}_{}".format(inst_name, constr_name),
                    main=new_model.rootAssembly.surfaces[
                        "{}_{}".format(inst_name, constr_obj.main[0])
                    ],
                    secondary=new_model.rootAssembly.surfaces[
                        "{}_{}".format(inst_name, constr_obj.secondary[0])
                    ],
                    constraintEnforcement=SURFACE_TO_SURFACE,
                    tieRotations=ON,
                    positionToleranceMethod=COMPUTED,
                    thickness=ON,
                    adjust=OFF,
                )


def copy_connector_assignments_from_original_model(
    src_model, new_model, num_copies, part_name
):
    """Copy the connector assignments from the original model to the new model.

    Args:
        src_model: Model object.
        new_model: Model object.
        num_copies: Number of copies.
        part_name: Name of the part.
    """
    num_original_model_wires = len(src_model.rootAssembly.sectionAssignments)

    for i in range(num_copies):
        inst_name = "{}-{}".format(part_name, i + 1)
        instance_nodes = new_model.rootAssembly.instances[inst_name].nodes
        for j, _ in enumerate(src_model.rootAssembly.sectionAssignments):
            constraints_to_consider = MAP_COUPLINGS[j + 1]
            constraint_names = tuple(
                [
                    "{}_{}".format(inst_name, constraints_to_consider[k])
                    for k in range(len(constraints_to_consider))
                ]
            )

            nodes = []
            for constraint_name in constraint_names:
                control_point_nodes = new_model.rootAssembly.allInternalSets[
                    new_model.constraints[constraint_name].controlPoint[0]
                ].nodes
                try:
                    nodes.append(instance_nodes[control_point_nodes[0].label - 1])
                except:
                    closest_node = instance_nodes.getClosest(
                        coordinates=control_point_nodes[0].coordinates,
                        numToFind=2,
                        searchTolerance=0.01,
                    )
                    nodes.append(closest_node[0] if j % 2 != 0 else closest_node[1])
            nodes = [nodes[0], None] if len(nodes) == 1 else nodes
            new_model.rootAssembly.WirePolyLine(
                points=((nodes[0], nodes[1]),), mergeType=IMPRINT, meshable=False
            )

            new_model.rootAssembly.DatumCsysByThreePoints(
                name="{}_{}".format(
                    inst_name,
                    src_model.rootAssembly.features[
                        src_model.rootAssembly.features.keys()[j + 1]
                    ].name,
                ),
                coordSysType=CARTESIAN,
                line1=(1.0, 0.0, 0.0),
                line2=(0.0, 1.0, 0.0),
                origin=(0.0, 0.0, 0.0),
            )

            index = new_model.rootAssembly.edges[0].index
            edge_array = new_model.rootAssembly.edges[index : index + 1]
            region = rts.Region(edges=edge_array)

            datum_id = new_model.rootAssembly.datums.keys()[-1]

            edge_feature_str = new_model.rootAssembly.edges[index].featureName
            edge_feature_num = int(edge_feature_str.split("-")[1])
            section_assignment_number = (
                edge_feature_num - 1
            ) % num_original_model_wires

            section_assignment_name = src_model.rootAssembly.sectionAssignments[
                section_assignment_number
            ].sectionName

            new_model.rootAssembly.SectionAssignment(
                region=region, sectionName=section_assignment_name
            )

            new_model.rootAssembly.ConnectorOrientation(
                region=new_model.rootAssembly.allInternalSets[
                    new_model.rootAssembly.allInternalSets.keys()[-1]
                ],
                axis1=AXIS_1,
                axis2=AXIS_1,
                angle1=0.0,
                angle2=0.0,
                localCsys1=new_model.rootAssembly.datums[datum_id],
                orient2sameAs1=True,
            )


def create_tying_constraints_between_instances(new_model, num_copies, part_name):
    """Create tie constraints between the boundary conditions surfaces (1 attached to 4, 2 attached to 3).

    Args:
        new_model: Model object.
        num_copies: Number of copies.
        part_name: Name of the part.
    """
    # add tie constraints between the boundary conditions surfaces (1 attached to 4, 2 attached to 3)
    for i in range(num_copies - 1):
        main_inst_name = "{}-{}".format(part_name, i + 1)
        secondary_inst_name = "{}-{}".format(part_name, i + 2)
        new_model.Tie(
            name="{}-TIE-{}-1".format(main_inst_name, secondary_inst_name),
            main=new_model.rootAssembly.surfaces["{}_Disp-BC-4".format(main_inst_name)],
            secondary=new_model.rootAssembly.surfaces[
                "{}_Disp-BC-1".format(secondary_inst_name)
            ],
            constraintEnforcement=SURFACE_TO_SURFACE,
            tieRotations=ON,
            positionToleranceMethod=COMPUTED,
            thickness=ON,
            adjust=OFF,
        )
        new_model.Tie(
            name="{}-TIE-{}-2".format(main_inst_name, secondary_inst_name),
            main=new_model.rootAssembly.surfaces["{}_Disp-BC-3".format(main_inst_name)],
            secondary=new_model.rootAssembly.surfaces[
                "{}_Disp-BC-2".format(secondary_inst_name)
            ],
            constraintEnforcement=SURFACE_TO_SURFACE,
            tieRotations=ON,
            positionToleranceMethod=COMPUTED,
            thickness=ON,
            adjust=OFF,
        )


def define_frequency_extraction_job(new_model_name, num_eigenmodes=50, solve=True):
    """Define the frequency extraction job.

    Args:
        new_model_name: Name of the new model.
        num_eigenvalues: Number of eigenvalues to extract.
        solve: Whether to solve the job.
    """
    mdb.models[new_model_name].FrequencyStep(
        name="Frequency Step1",
        previous="Initial",
        numEigen=num_eigenmodes,
        eigensolver=LANCZOS,
        normalization=MASS,
        acousticCoupling=AC_PROJECTION,
    )
    mdb.Job(name="FrequencyExtraction", model=new_model_name, resultsFormat=BOTH)
    mdb.jobs["FrequencyExtraction"].setValues(
        multiprocessingMode=MPI, numCpus=4, numDomains=4, numThreadsPerMpiProcess=1
    )

    start_time = time.time()
    if solve:
        mdb.jobs["FrequencyExtraction"].submit()
        mdb.jobs["FrequencyExtraction"].waitForCompletion()

    print("  - Frequency extraction took: {:.2f}s".format(time.time() - start_time))


def create_rail_top_union_set(new_model, num_copies, part_name="PART-1"):
    """Create a union set of the rail top sets.

    Args:
        new_model: Model object.
        num_copies: Number of copies.
        part_name: Name of the part.
    """
    new_asm = new_model.rootAssembly
    for set_name in ["RailTop-1 (847)-1", "RailTop-2 (849)-1"]:
        sets = []
        for i in range(num_copies):
            inst_name = "{}-{}".format(part_name, i + 1)
            sets.append(new_asm.sets["{}_{}".format(inst_name, set_name)])
        new_asm.SetByBoolean(name=set_name, sets=sets, operation=UNION)
        for set in sets:
            del set


def copy_part_with_assembly_sets_surfaces(
    path,
    source_model_name="Linear_Flextrack",
    part_name="PART-1",
    new_model_name="Linear_Flextrack_pattern",
    num_copies=5,
    solve=True,
):
    """Copy the part with assembly sets and surfaces.

    Args:
        path: Path to the source model.
        source_model_name: Name of the source model.
        part_name: Name of the part.
        new_model_name: Name of the new model.
        num_copies: Number of copies.
        solve: Whether to submit the job and do the analysis.
    """
    # Open source model
    start_time = time.time()
    src_model = open_source_model(path, source_model_name)
    src_part = src_model.parts[part_name]
    src_asm = src_model.rootAssembly
    x_coords = [node.coordinates[0] for node in src_part.nodes]
    track_length = max(x_coords) - min(x_coords)

    add_tying_surfaces_for_replication(src_model, part_name, num_copies)

    new_model = define_new_model(new_model_name, part_name, src_part)
    part_copy = new_model.parts[part_name]
    new_asm = new_model.rootAssembly

    define_model_materials(new_model)
    copy_part_and_translate(new_asm, part_name, part_copy, track_length, num_copies)
    copy_sections_from_original_model(src_model, new_model, num_copies, part_name)
    copy_sets_and_internal_sets_from_original_model(
        src_asm, new_asm, num_copies, part_name
    )
    copy_surfaces_from_original_model(
        src_asm, new_asm, src_model, num_copies, part_name
    )
    define_and_remove_redundant_bcs(src_model, new_model, num_copies, part_name)
    copy_constraints_from_original_model(src_model, new_model, num_copies, part_name)
    copy_connector_assignments_from_original_model(
        src_model, new_model, num_copies, part_name
    )
    create_tying_constraints_between_instances(new_model, num_copies, part_name)
    create_rail_top_union_set(new_model, num_copies, part_name)
    define_frequency_extraction_job(new_model_name, solve=solve)
    mdb.saveAs(path)


if __name__ == "__main__":
    start_time = time.time()
    copy_part_with_assembly_sets_surfaces(
        path="data/simpack/fe-data/Linear_FlexTrack_2023.cae",
        source_model_name="Linear_Flextrack",
        part_name="PART-1",
        new_model_name="Linear_Flextrack_pattern",
        num_copies=1,
        solve=False,
    )
    print("  - Total time: {:.2f}s".format(time.time() - start_time))

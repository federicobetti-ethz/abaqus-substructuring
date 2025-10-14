from abaqus import mdb, session
from abaqusConstants import *
import regionToolset
import mesh
import time

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse import save_npz

model_name = 'BlockModel'
part_name = 'Block'
material_name = 'Steel'
section_name = 'BlockSection'
assembly_name = 'Assembly'
instance_name = part_name + '-1'

L = 10.0
H = 0.2
W = 0.2

num_substructures = 4

seed_size = 0.05

E = 210e9
nu = 0.3
rho = 7850

cae_file = 'block_model.cae'
input_file = 'block_model.inp'


if model_name in mdb.models:
    del mdb.models[model_name]

model = mdb.Model(name=model_name)
s = model.ConstrainedSketch(name='__profile__', sheetSize=2.0)
s.rectangle(point1=(0.0, 0.0), point2=(L, H))

part = model.Part(name=part_name, dimensionality=THREE_D, type=DEFORMABLE_BODY)
part.BaseSolidExtrude(sketch=s, depth=W)

material = model.Material(name=material_name)
material.Elastic(table=((E, nu),))
material.Density(table=((rho,),))

model.HomogeneousSolidSection(
    name=section_name,
    material=material_name,
    thickness=None
)

cells = part.cells[:]
region = regionToolset.Region(cells=cells)
part.SectionAssignment(region=region, sectionName=section_name)

assembly = model.rootAssembly
instance = assembly.Instance(name=instance_name, part=part, dependent=ON)

faces = instance.faces
leftFaces = faces.getByBoundingBox(xMin=-1e-6, xMax=1e-6)

if len(leftFaces) == 0:
    leftFaces = [f for f in faces if abs(f.pointOn[0]) < 1e-6]

assembly.Set(name='LeftFace', faces=leftFaces)
part.seedPart(size=seed_size, deviationFactor=0.1, minSizeFactor=0.1)

elemType1 = mesh.ElemType(elemCode=C3D8, elemLibrary=STANDARD)
part.setElementType(regions=(part.cells[:],), elemTypes=(elemType1,))
part.generateMesh()

assembly.regenerate()

leftSet = assembly.sets['LeftFace']
bc_labels = [node.label for node in leftSet.nodes]
model.EncastreBC(name='BC_Clamped_Left', createStepName='Initial', region=leftSet)

for i in range(num_substructures):
    internal_nodes = []
    surface_nodes = []
    for node in instance.nodes:
        if node.coordinates[0] > i * L / num_substructures and node.coordinates[0] < (i + 1) * L / num_substructures:
            internal_nodes.append(node.label)
        elif np.isclose(node.coordinates[0], i * L / num_substructures) or np.isclose(node.coordinates[0], (i + 1) * L / num_substructures):
            surface_nodes.append(node.label)
    internal_node_set = instance.nodes.sequenceFromLabels(labels=internal_nodes)
    surface_node_set = instance.nodes.sequenceFromLabels(labels=surface_nodes)
    assembly.Set(name=f'Substructure{i}InternalNodes', nodes=internal_node_set)
    assembly.Set(name=f'Substructure{i}SurfaceNodes', nodes=surface_node_set)

model.FrequencyStep(name='ExtractMatrices', previous='Initial', numEigen=1)

mdb.saveAs(pathName=cae_file)
mdb.Job(name='BlockMatrixJob', model=model_name).writeInput(consistencyChecking=OFF)

input_file = 'BlockMatrixJob.inp'
with open(input_file, 'r') as f:
    lines = f.readlines()

target_line = "*Step, name=ExtractMatrices, nlgeom=NO, perturbation"
for i, line in enumerate(lines):
    if target_line in line:
        lines[i] = "*Step, name=ExtractMatrices\n"

target_line = '*Frequency, eigensolver=Lanczos, sim, acoustic coupling=on, normalization=mass'
for i, line in enumerate(lines):
    if target_line in line:
        lines[i] = "*MATRIX GENERATE, STIFFNESS, MASS\n"

target_line = "1, , , , ,"
for i, line in enumerate(lines):
    if target_line in line:
        lines[i] = "*MATRIX OUTPUT, STIFFNESS, MASS, FORMAT=MATRIX INPUT\n"

end_line = None
target_line = "** OUTPUT REQUESTS"
for i, line in enumerate(lines):
    if target_line in line:
        end_line = i-1

num_lines_to_erase = len(lines) - 1 - end_line
for idx in range(num_lines_to_erase):
    lines.pop(-2)

with open(input_file, 'w') as f:
    f.writelines(lines)

time.sleep(5)

mdb.JobFromInputFile(name='BlockMatrixJob', inputFileName=input_file)
mdb.jobs['BlockMatrixJob'].submit()
mdb.jobs['BlockMatrixJob'].waitForCompletion()

bc_dofs = []
n_dofs_per_node = 3
for label in bc_labels:
    label = label - 1
    for dof in range(n_dofs_per_node):
        bc_dofs.append(label * n_dofs_per_node + dof)

time.sleep(5)

data = np.loadtxt('BlockMatrixJob_STIF1.mtx', delimiter=',')

node_i = data[:, 0].astype(int) - 1
dof_i  = data[:, 1].astype(int) - 1
node_j = data[:, 2].astype(int) - 1
dof_j  = data[:, 3].astype(int) - 1
vals   = data[:, 4]

rows = node_i * n_dofs_per_node + dof_i
cols = node_j * n_dofs_per_node + dof_j

K = coo_matrix((vals, (rows, cols))).tocsr()
K = 0.5 * (K + K.T)

data = np.loadtxt('BlockMatrixJob_MASS1.mtx', delimiter=',')

node_i = data[:, 0].astype(int) - 1
dof_i  = data[:, 1].astype(int) - 1
node_j = data[:, 2].astype(int) - 1
dof_j  = data[:, 3].astype(int) - 1
vals   = data[:, 4]

n_dofs_per_node = 3
rows = node_i * n_dofs_per_node + dof_i
cols = node_j * n_dofs_per_node + dof_j

M = coo_matrix((vals, (rows, cols))).tocsr()

for i in range(num_substructures):
    internal_nodes = assembly.sets[f'Substructure{i}InternalNodes']
    surface_nodes = assembly.sets[f'Substructure{i}SurfaceNodes']
    internal_node_labels = [node.label for node in internal_nodes.nodes]
    surface_node_labels = [node.label for node in surface_nodes.nodes]

    internal_node_labels = np.array(internal_node_labels) - 1
    surface_node_labels = np.array(surface_node_labels) - 1
    
    internal_dofs = internal_node_labels[:, np.newaxis] * n_dofs_per_node + np.arange(n_dofs_per_node)
    surface_dofs = surface_node_labels[:, np.newaxis] * n_dofs_per_node + np.arange(n_dofs_per_node)
    
    internal_dofs = internal_dofs.flatten()
    surface_dofs = surface_dofs.flatten()

    dofs_to_eliminate = np.intersect1d(surface_dofs, bc_dofs)
    not_bc_dofs = np.setdiff1d(np.arange(M.shape[0]), dofs_to_eliminate)
    surface_dofs = np.intersect1d(surface_dofs, not_bc_dofs)

    Mr_ii = M[np.ix_(internal_dofs, internal_dofs)]
    Kr_ii = K[np.ix_(internal_dofs, internal_dofs)]
    Mr_bb = M[np.ix_(surface_dofs, surface_dofs)]
    Kr_bb = K[np.ix_(surface_dofs, surface_dofs)]

    Kr_ib = K[np.ix_(internal_dofs, surface_dofs)]
    Kr_bi = K[np.ix_(surface_dofs, internal_dofs)]

    save_npz(f'Mr_{i}_ii.npz', Mr_ii)
    save_npz(f'Kr_{i}_ii.npz', Kr_ii)
    save_npz(f'Mr_{i}_bb.npz', Mr_bb)
    save_npz(f'Kr_{i}_bb.npz', Kr_bb)
    save_npz(f'Kr_{i}_ib.npz', Kr_ib)
    save_npz(f'Kr_{i}_bi.npz', Kr_bi)

    np.savetxt(f'internal_dofs_{i}.txt', internal_dofs)
    np.savetxt(f'surface_dofs_{i}.txt', surface_dofs)

with open('bc_dofs.txt', 'w') as f:
    for dof in bc_dofs:
        f.write(f"{dof}\n")

with open('num_total_dofs.txt', 'w') as f:
    f.write(f"{M.shape[0]}\n")

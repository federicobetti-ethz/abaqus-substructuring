from abaqus import mdb, session
from abaqusConstants import *
import regionToolset
import mesh
import time

import numpy as np
from scipy.io import mmread
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import eigsh

model_name = 'BlockModel'
part_name = 'Block'
material_name = 'Steel'
section_name = 'BlockSection'
assembly_name = 'Assembly'
instance_name = part_name + '-1'

L = 1.0
H = 0.2
W = 0.2

seed_size = 0.1

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
assembly.Set(name='WholeInstance', nodes=instance.nodes)

part.seedPart(size=seed_size, deviationFactor=0.1, minSizeFactor=0.1)

elemType1 = mesh.ElemType(elemCode=C3D8, elemLibrary=STANDARD)
part.setElementType(regions=(part.cells[:],), elemTypes=(elemType1,))
part.generateMesh()

assembly.regenerate()

leftSet = assembly.sets['LeftFace']
bc_labels = [node.label for node in leftSet.nodes]
np.savetxt('leftSet.txt', bc_labels)
# model.EncastreBC(name='BC_Clamped', createStepName='Initial', region=leftSet)

model.FrequencyStep(name='ExtractMatrices', previous='Initial', numEigen=1)

mdb.saveAs(pathName=cae_file)
mdb.Job(name='BlockMatrixJob', model=model_name).writeInput(consistencyChecking=OFF)

input_file = 'BlockMatrixJob.inp'
with open(input_file, 'r') as f:
    lines = f.readlines()

target_line = '*Frequency, eigensolver=Lanczos, sim, acoustic coupling=on, normalization=mass'
for i, line in enumerate(lines):
    if target_line in line:
        lines[i] = "*MATRIX GENERATE, STIFFNESS, MASS\n"

target_line = "1, , , , ,"
for i, line in enumerate(lines):
    if target_line in line:
        lines[i] = "*MATRIX OUTPUT, STIFFNESS, MASS, FORMAT=MATRIX INPUT\n"

with open(input_file, 'w') as f:
    f.writelines(lines)

time.sleep(5)

mdb.JobFromInputFile(name='BlockMatrixJob', inputFileName=input_file)
mdb.jobs['BlockMatrixJob'].submit()
mdb.jobs['BlockMatrixJob'].waitForCompletion()

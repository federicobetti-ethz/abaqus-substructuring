print("Computed modes")

import numpy as np
from odbAccess import Odb
from abaqusConstants import *
import regionToolset

E = 210e9
nu = 0.3
rho = 7850

odb = session.Odb(name='BlockMatrixJob', path='BlockMatrixJob_from_python.odb')

mdb = openMdb(pathName='block_model.cae')
model = mdb.models['BlockModel']

scat = odb.SectionCategory(name="solid", description="Test")

nodes = model.rootAssembly.instances['Block-1'].nodes
elements = model.rootAssembly.instances['Block-1'].elements
material_name = 'Steel'
material = model.materials[material_name]

odb_material = odb.Material(name=material_name)
odb_material.Elastic(table=((E, nu),))
odb_material.Density(table=((rho,),))

section_name = 'SteelSection'
section = odb.HomogeneousSolidSection(
    name=section_name,
    material=material_name,
    thickness=2.0
)

part = odb.Part(name='PART-1', embeddedSpace=THREE_D, type=DEFORMABLE_BODY)
nodeData = []
for node in nodes:
    nodeData.append((node.label, *node.coordinates))
part.addNodes(nodeData=nodeData, nodeSetName='NodeSet-1')

elementData = []
for element in elements:
    element_nodes = element.getNodes()
    element_nodes_labels = [node.label for node in element_nodes]
    elementData.append((element.label, *element_nodes_labels))
part.addElements(
    elementData=elementData,
    type="C3D8",
    elementSetName='ElementSet-1',
    sectionCategory=scat
)

instance = odb.rootAssembly.Instance(name='PART-1-1', object=part)
element_set = instance.ElementSetFromElementLabels(
    name="AllElements-1",
    elementLabels=[elementData[0] for elementData in elementData]
)
instance.assignSection(region=element_set, section=section)
instance.NodeSetFromNodeLabels(
    name="AllNodes-1",
    nodeLabels=[nodeData[0] for nodeData in nodeData]
)

step = odb.Step(
    name='EigenStep',
    description='Eigenvalue Modes',
    domain=MODAL,
)

nodeLabels = [nodeDataPoint[0] for nodeDataPoint in nodeData]

num_substructures = 4
for i in range(num_substructures):
    vecs = np.loadtxt(f'eigenmodes_{i}.txt')
    freqs = np.loadtxt(f'eigenfrequencies_{i}.txt')

    nDofs, nModes = vecs.shape

    for j in range(nModes):
        frame = step.Frame(mode=j, frequency=0.1)
        disp = frame.FieldOutput(
            name='U',
            description=f'Mode {j+1}',
            type=VECTOR,
            validInvariants=(MAGNITUDE,)
        )
        labels = tuple(nodeLabels)
        data = vecs[:, j]
        data = data.reshape(-1, 3)
        data = tuple(tuple(data[k]) for k in range(data.shape[0]))
        disp.addData(position=NODAL, instance=instance, labels=labels, data=data)

odb.save()
odb.close()

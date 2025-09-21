import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import eigsh

K = np.loadtxt("BlockMatrixJob_STIF1.mtx", delimiter=',', usecols=(0, 2, 4))
M = np.loadtxt("BlockMatrixJob_MASS1.mtx", delimiter=',', usecols=(0, 2, 4))

bc_labels = np.loadtxt('leftSet.txt')
bc_labels = bc_labels.astype(int)-1
bc_dofs = []
n_dofs_per_node = 3
for label in bc_labels:
    for dof in range(n_dofs_per_node):
        bc_dofs.append(label * n_dofs_per_node + dof)

data = np.loadtxt('BlockMatrixJob_STIF1.mtx', delimiter=',')

node_i = data[:,0].astype(int) - 1
dof_i  = data[:,1].astype(int) - 1
node_j = data[:,2].astype(int) - 1
dof_j  = data[:,3].astype(int) - 1
vals   = data[:,4]

n_dofs_per_node = 3
rows = node_i * n_dofs_per_node + dof_i
cols = node_j * n_dofs_per_node + dof_j

K = coo_matrix((vals, (rows, cols))).tocsr()
K = K.todense()

K[np.ix_(bc_dofs, bc_dofs)] = 0
K[bc_dofs, bc_dofs] = 1

data = np.loadtxt('BlockMatrixJob_MASS1.mtx', delimiter=',')

node_i = data[:,0].astype(int) - 1
dof_i  = data[:,1].astype(int) - 1
node_j = data[:,2].astype(int) - 1
dof_j  = data[:,3].astype(int) - 1
vals   = data[:,4]

n_dofs_per_node = 3
rows = node_i * n_dofs_per_node + dof_i
cols = node_j * n_dofs_per_node + dof_j

M = coo_matrix((vals, (rows, cols))).tocsr()
M = M.todense()

vals, vecs = eigsh(K, 20, M=M)
freqs = np.sqrt(vals) / (2 * np.pi)

np.savetxt('eigenmodes.txt', vecs)
np.savetxt('eigenfreqs.txt', freqs)
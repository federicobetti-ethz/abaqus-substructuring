"""Main script to run the substructuring method."""

import numpy as np
from scipy.sparse import load_npz
from scipy.sparse.linalg import eigsh, spsolve
import subprocess
import os

cmd = f'abaqus cae noGUI={os.path.join(os.path.dirname(__file__), "example.py")}'
subprocess.run(cmd, shell=True, check=True)

num_total_dofs = np.loadtxt("num_total_dofs.txt").astype(int)

num_substructures = 4
for i in range(num_substructures):
    Kr_ii = load_npz(f"Kr_{i}_ii.npz")
    Mr_ii = load_npz(f"Mr_{i}_ii.npz")
    Kr_bb = load_npz(f"Kr_{i}_bb.npz")
    Mr_bb = load_npz(f"Mr_{i}_bb.npz")
    Kr_ib = load_npz(f"Kr_{i}_ib.npz")
    Kr_bi = load_npz(f"Kr_{i}_bi.npz")

    internal_dofs = np.loadtxt(f"internal_dofs_{i}.txt").astype(int)
    surface_dofs = np.loadtxt(f"surface_dofs_{i}.txt").astype(int)

    eigvals, eigvecs = eigsh(Kr_ii, 20, M=Mr_ii, which="SM")
    static_constraint_modes = -spsolve(Kr_ii.tocsc(), Kr_ib.tocsc()).toarray()
    eigvecs = np.concatenate((eigvecs, static_constraint_modes), axis=1)
    eigenmodes = np.zeros((num_total_dofs, eigvecs.shape[1]))
    eigenmodes[internal_dofs, :] = eigvecs
    freqs = np.sqrt(eigvals) / (2 * np.pi)

    np.savetxt(f"eigenmodes_{i}.txt", eigenmodes)
    np.savetxt(f"eigenfrequencies_{i}.txt", freqs)

cmd = f'abaqus cae noGUI={os.path.join(os.path.dirname(__file__), "odb_example.py")}'
subprocess.run(cmd, shell=True, check=True)

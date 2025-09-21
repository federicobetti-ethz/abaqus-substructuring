# abaqus-substructuring
This is a small repository to try to do the following:

1. create an Abaqus model with scripting
2. extract stiffness and mass matrices
3. do modal analysis with a custom solver in Python
4. write results to an output database
5. visualize results in Abaqus CAE and verify they are correct

One should execute the following commands to test the current scripts:

```
abaqus cae noGUI=example.py
python eigenmodes_example.py
abaqus cae noGUI=odb_example.py
```
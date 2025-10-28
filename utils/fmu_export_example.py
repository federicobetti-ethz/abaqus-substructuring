"""Example of how to export an FMU from the simulator."""

from fmpy import dump, simulate_fmu
import matplotlib.pyplot as plt
import time
import numpy as np

fmu = "out/results_0/models/Train.output/Train.fmu"
dump(fmu)

start = time.time()
result = simulate_fmu(fmu, stop_time=100.0)
print(f"Simulation finished in {time.time() - start} seconds")

times = np.array([result[i][0] for i in range(1, len(result))])
positions = np.array([result[i][1] for i in range(1, len(result))])
plt.grid(True)
plt.xlabel("Time (s)")
plt.ylabel("Position (m)")
plt.title("Train Position vs Time")
plt.plot(times, positions)
plt.show()

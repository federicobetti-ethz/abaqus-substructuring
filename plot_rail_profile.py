import matplotlib.pyplot as plt
import numpy as np

# Rail profile points from rail_mesh_refinement.py (reduced bottom line)
rail_profile_points = [
    (1.316689, 423.774E-03),
    (1.307883, 418.203E-03),
    (1.297689, 414.324E-03),
    (1.293897, 410.876E-03),
    (1.291245, 403.458E-03),
    (1.288447, 385.816E-03),
    (1.288039, 375.928E-03),
    (1.288039, 365.147E-03),
    (1.288039, 343.928E-03),
    (1.28855, 332.865E-03),
    (1.290056, 322.022E-03),
    (1.292114, 313.325E-03),
    (1.29636, 308.967E-03),
    (1.314243, 302.464E-03),
    (1.334994, 296.545E-03),
    (1.351075, 295.394E-03),
    (1.354789, 291.404E-03),
    (1.354789, 285.628E-03),
    (1.340813, 283.628E-03),
    (1.279789, 283.628E-03),
    (1.217987, 283.628E-03),
    (1.206789, 283.628E-03),
    (1.204789, 291.404E-03),
    (1.208504, 295.394E-03),
    (1.217572, 296.043E-03),
    (1.224584, 296.545E-03),
    (1.244588, 302.192E-03),
    (1.263218, 308.967E-03),
    (1.267464, 313.325E-03),
    (1.269523, 322.022E-03),
    (1.271028, 332.863E-03),
    (1.271539, 343.928E-03),
    (1.271539, 365.147E-03),
    (1.271539, 375.928E-03),
    (1.271131, 385.814E-03),
    (1.268339, 403.458E-03),
    (1.265682, 410.873E-03),
    (1.26189, 414.324E-03),
    (1.25177, 418.175E-03),
    (1.244818, 420.82E-03),
    (1.243191, 429.775E-03),
    (1.243773, 441.367E-03),
    (1.246324, 448.472E-03),
    (1.253665, 453.342E-03),
    (1.269539, 455.453E-03),
    (1.279789, 455.628E-03),
    (1.290039, 455.453E-03),
    (1.305914, 453.342E-03),
    (1.312847, 448.992E-03),
    (1.315806, 441.367E-03),
    (1.316355, 430.434E-03),
]


x_coords = [point[0] for point in rail_profile_points]
y_coords = [point[1] for point in rail_profile_points]

plt.figure(figsize=(15, 10))

plt.plot(x_coords + [x_coords[0]], y_coords + [y_coords[0]], 'b-', linewidth=2, label='Rail Profile')
plt.scatter(x_coords, y_coords, color='red', s=100, zorder=5, edgecolors='black', linewidth=2)

for i, (x, y) in enumerate(rail_profile_points):
    plt.annotate(f'{i}', 
                (x, y), 
                xytext=(0, 0), 
                textcoords='offset points',
                fontsize=8,
                ha='center', va='center',
                color='white', fontweight='bold')

plt.grid(True, alpha=0.3)
plt.xlabel('X-coordinate (m)', fontsize=12)
plt.ylabel('Y-coordinate (m)', fontsize=12)
plt.title(f'Rail Profile\n({len(rail_profile_points)} points total)', fontsize=14, fontweight='bold')

plt.axis('equal')

plt.xlim(min(x_coords) - 0.05, max(x_coords) + 0.05)
plt.ylim(min(y_coords) - 0.02, max(y_coords) + 0.02)
plt.legend(loc='upper right')

textstr = f'Rail Profile Analysis:\nTotal points: {len(rail_profile_points)}\nX range: {min(x_coords):.3f} to {max(x_coords):.3f} m\nY range: {min(y_coords):.3f} to {max(y_coords):.3f} m\n\nCentral refinement region:\nZ = {0.205:.3f} to {0.41:.3f} m'
props = dict(boxstyle='round', facecolor='wheat', alpha=0.8)
plt.text(0.02, 0.98, textstr, transform=plt.gca().transAxes, fontsize=10,
         verticalalignment='top', bbox=props)

plt.tight_layout()
plt.show()
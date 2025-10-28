import matplotlib.pyplot as plt
import numpy as np

left_half_points = [
    (1.263218, 283.628e-03),
    (1.244588, 283.628e-03),
    (1.224584, 283.628e-03),
    (1.204789, 285.628e-03),
    (1.204789, 291.404e-03),
    (1.208504, 295.394e-03),
    (1.224584, 296.545e-03),
    (1.244588, 302.192e-03),
    (1.263218, 308.967e-03),
    (1.267464, 313.325e-03),
    (1.269523, 322.022e-03),
    (1.271539, 343.928e-03),
    (1.271539, 359.538e-03),
    (1.271131, 385.814e-03),
    (1.268339, 403.458e-03),
    (1.265682, 410.873e-03),
    (1.26189, 414.324e-03),
    (1.25177, 418.175e-03),
    (1.244818, 420.82e-03),
    (1.243191, 429.775e-03),
    (1.243773, 441.367e-03),
    (1.246324, 448.472e-03),
    (1.253665, 453.342e-03),
    (1.269539, 455.453e-03),
]

center_top = [(1.279789, 455.628e-03)]
center_bottom = [(1.279789, 283.628e-03)]


def reflect_point(point):
    x, y = point
    center_x = center_top[0][0]
    new_x = 2 * center_x - x
    return (new_x, y)


right_half_points = [reflect_point(p) for p in reversed(left_half_points)]
rail_profile_points = left_half_points + center_top + right_half_points + center_bottom


x_coords = [point[0] for point in rail_profile_points]
y_coords = [point[1] for point in rail_profile_points]

plt.figure(figsize=(15, 10))

plt.plot(
    x_coords + [x_coords[0]],
    y_coords + [y_coords[0]],
    "b-",
    linewidth=2,
    label="Rail Profile",
)
plt.scatter(
    x_coords, y_coords, color="red", s=100, zorder=5, edgecolors="black", linewidth=2
)

for i, (x, y) in enumerate(rail_profile_points):
    plt.annotate(
        f"{i}",
        (x, y),
        xytext=(0, 0),
        textcoords="offset points",
        fontsize=8,
        ha="center",
        va="center",
        color="white",
        fontweight="bold",
    )

plt.grid(True, alpha=0.3)
plt.xlabel("X-coordinate (m)", fontsize=12)
plt.ylabel("Y-coordinate (m)", fontsize=12)
plt.title(
    f"Rail Profile\n({len(rail_profile_points)} points total)",
    fontsize=14,
    fontweight="bold",
)

plt.axis("equal")

plt.xlim(min(x_coords) - 0.05, max(x_coords) + 0.05)
plt.ylim(min(y_coords) - 0.02, max(y_coords) + 0.02)
plt.legend(loc="upper right")

textstr = f"Rail Profile Analysis:\nTotal points: {len(rail_profile_points)}\nX range: {min(x_coords):.3f} to {max(x_coords):.3f} m\nY range: {min(y_coords):.3f} to {max(y_coords):.3f} m\n\nCentral refinement region:\nZ = {0.205:.3f} to {0.41:.3f} m"
props = dict(boxstyle="round", facecolor="wheat", alpha=0.8)
plt.text(
    0.02,
    0.98,
    textstr,
    transform=plt.gca().transAxes,
    fontsize=10,
    verticalalignment="top",
    bbox=props,
)

plt.tight_layout()
plt.show()

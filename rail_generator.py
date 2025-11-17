import gmsh
import numpy as np


def generate_rail_sweep(
    txt_file: str,
    rail_profile_points: list,
    segment_length: float = 100.0,
    overlap: float = 12.0,
    min_segment_length: float = 20.0,
    step_size: int = 1,
    step_file: str = "rail_sweep.step",
    model_name: str = "RailSweep",
):
    """
    Generate a rail by sweeping a profile along a spline, split into multiple segments based on arc length.

    Args:
        txt_file: Path to trajectory data file
        rail_profile_points: List of (y, z) tuples defining the rail profile
        segment_length: Length of each segment along the curve in meters (default: 100.0)
        overlap: Overlap between adjacent segments in meters (default: 12.0)
        min_segment_length: Minimum segment length to keep, otherwise skip last segment (default: 20.0)
        step_size: Downsampling factor for spline points (default: 10)
                   Lower values (2-3) give smoother curves but more data.
                   Use step_size=2 or 3 for smooth results.
        step_file: Output STEP file path
        model_name: Name for the GMSH model
    """
    gmsh.initialize()
    gmsh.model.add(model_name)

    y_coords = [y for y, _ in rail_profile_points]
    z_coords = [z for _, z in rail_profile_points]
    y_center = (min(y_coords) + max(y_coords)) / 2
    z_center = (min(z_coords) + max(z_coords)) / 2
    rail_profile_points = [(y - y_center, z - z_center) for y, z in rail_profile_points]

    ss, xs, ys, zs = [], [], [], []
    with open(txt_file, "r") as f:
        lines = f.readlines()

        for idx, line in enumerate(lines):
            if line.startswith('"x(s)"'):
                xline = idx
            elif line.startswith('"y(s)"'):
                yline = idx
            elif line.startswith('"z(s)"'):
                zline = idx
            elif line.startswith('"Superelevation u(s)"'):
                superelevation_line = idx

        for line in lines[xline + 3 : yline]:
            parts = line.split(",")
            ss.append(float(parts[0]))
            xs.append(float(parts[1]))

        for line in lines[yline + 3 : zline]:
            parts = line.split(",")
            ys.append(float(parts[1]))

        for line in lines[zline + 3 : superelevation_line]:
            parts = line.split(",")
            zs.append(float(parts[1]))

    min_len = min(len(ss), len(xs), len(ys), len(zs))
    ss = np.array(ss[:min_len])
    xs = np.array(xs[:min_len])
    ys = np.array(ys[:min_len])
    zs = np.array(zs[:min_len])

    if step_size > 1:
        indices_ds = np.arange(0, len(xs), step_size)
        ss_ds = ss[indices_ds]
        xs_ds = xs[indices_ds]
        ys_ds = ys[indices_ds]
        zs_ds = zs[indices_ds]
    else:
        ss_ds = ss
        xs_ds = xs
        ys_ds = ys
        zs_ds = zs

    total_length = ss[-1] - ss[0]
    print(f"Total track length: {total_length:.2f} m")

    if total_length <= segment_length:
        print(
            f"Track shorter than segment length ({total_length:.2f} m < {segment_length:.2f} m), creating single segment"
        )
        start_idx = 0
        end_idx = len(ss_ds)
        segments = [(start_idx, end_idx, ss[0], ss[-1])]
    else:
        net_segment_length = segment_length - overlap
        segments = []

        current_s = ss[0]
        seg_idx = 0

        while current_s < ss[-1]:
            segment_start_s = current_s
            segment_end_s = min(current_s + segment_length, ss[-1])

            remaining_after_this = ss[-1] - segment_end_s
            if remaining_after_this < min_segment_length and remaining_after_this > 0:
                segment_end_s = ss[-1]
                print(
                    f"Extending segment to cover remaining {segment_end_s - segment_start_s:.2f} m (to avoid {remaining_after_this:.2f} m short segment)"
                )

            start_idx = np.searchsorted(ss_ds, segment_start_s)
            end_idx = np.searchsorted(ss_ds, segment_end_s, side="right")

            if end_idx - start_idx < 2:
                break

            segments.append((start_idx, end_idx, segment_start_s, segment_end_s))

            if segment_end_s >= ss[-1]:
                break

            current_s += net_segment_length
            seg_idx += 1

    print(f"Creating {len(segments)} segments")

    volumes = []
    segment_metadata = []

    for seg_idx, (start_idx, end_idx, seg_start_s, seg_end_s) in enumerate(segments):
        print(
            f"Segment {seg_idx + 1}: s = {seg_start_s:.2f} to {seg_end_s:.2f} m (length = {seg_end_s - seg_start_s:.2f} m)"
        )

        seg_xs = xs_ds[start_idx:end_idx]
        seg_ys = ys_ds[start_idx:end_idx]
        seg_zs = zs_ds[start_idx:end_idx]

        seg_start_x, seg_start_y, seg_start_z = seg_xs[0], seg_ys[0], seg_zs[0]

        profile_pts_tags = []
        for y, z in rail_profile_points:
            profile_pts_tags.append(
                gmsh.model.occ.addPoint(seg_start_x, seg_start_y - y, seg_start_z - z)
            )

        profile_lines = []
        for i in range(len(profile_pts_tags)):
            pstart = i
            pend = (i + 1) % len(profile_pts_tags)
            profile_lines.append(
                gmsh.model.occ.addLine(profile_pts_tags[pstart], profile_pts_tags[pend])
            )

        cl = gmsh.model.occ.addCurveLoop(profile_lines)
        profile_face = gmsh.model.occ.addPlaneSurface([cl])

        segment_pts = []
        for i in range(len(seg_xs)):
            segment_pts.append(gmsh.model.occ.addPoint(seg_xs[i], seg_ys[i], seg_zs[i]))

        segment_spline = gmsh.model.occ.addSpline(segment_pts)
        wire_tag = gmsh.model.occ.addWire([segment_spline])
        sweep = gmsh.model.occ.addPipe([(2, profile_face)], wire_tag)

        for dim, tag in sweep:
            if dim == 3:
                volumes.append(tag)

    gmsh.model.occ.synchronize()

    gmsh.write(step_file)
    gmsh.finalize()

    for seg_idx, (_, _, seg_start_s, seg_end_s) in enumerate(segments, start=1):
        segment_metadata.append(
            {
                "idx": seg_idx,
                "s_start": float(seg_start_s),
                "s_end": float(seg_end_s),
            }
        )

    return segment_metadata

import argparse
import sys
from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
import numpy as np

NUM_LEVELS = 20
CONTOUR_LINES = 10
FIGSIZE_WIDTH = 14.0


def read_state_cache(filename):
    with open(filename, "r", encoding="utf-8") as f:
        lines = f.readlines()

    nx = ny = None
    for line in lines[:3]:
        if line.startswith("# NX"):
            parts = line.split()
            nx = int(parts[2])
            ny = int(parts[4])
            break

    data = []
    for line in lines:
        if line.startswith("#"):
            continue
        parts = line.strip().split()
        if len(parts) < 8:
            continue
        xi, eta, x, y, w, psi, ux, uy = map(float, parts[:8])
        data.append((xi, eta, x, y, w, psi, ux, uy))

    data = np.asarray(data, dtype=float)
    if nx is None or ny is None:
        nx = int(np.sqrt(len(data)))
        ny = nx

    return data, nx, ny


def display_prop_name(prop_name):
    if prop_name == "psi":
        return "ψ"
    return prop_name


def metadata_text(reynolds_number, nx_meta, ny_meta, u0):
    parts = []
    if reynolds_number is not None:
        parts.append(f"Re={reynolds_number:g}")
    if nx_meta is not None and ny_meta is not None:
        parts.append(f"Grid={nx_meta}x{ny_meta}")
    if u0 is not None:
        parts.append(f"u0={u0:g}")
    return " | ".join(parts)


def save_path(prop_name, reynolds_number, nx_meta, ny_meta, u0):
    images_dir = Path("Images")
    images_dir.mkdir(exist_ok=True)

    suffix_parts = []
    if reynolds_number is not None:
        suffix_parts.append(f"Re{reynolds_number:g}")
    if nx_meta is not None and ny_meta is not None:
        suffix_parts.append(f"N{nx_meta}x{ny_meta}")
    if u0 is not None:
        suffix_parts.append(f"u0{u0:g}")

    suffix = "_" + "_".join(suffix_parts) if suffix_parts else ""
    return images_dir / f"state_cache_{prop_name}{suffix}.png"


def plot_property(
    filename,
    prop_name,
    num_levels=None,
    contour_lines=None,
    figsize_width=None,
    plot_areas=True,
    plot_lines=True,
    linewidth=0.5,
    alpha=0.3,
    save_images=False,
    reynolds_number=None,
    nx_meta=None,
    ny_meta=None,
    u0=None,
):
    if num_levels is None:
        num_levels = NUM_LEVELS
    if contour_lines is None:
        contour_lines = CONTOUR_LINES
    if figsize_width is None:
        figsize_width = FIGSIZE_WIDTH

    figsize = (figsize_width, 6)
    data, nx, ny = read_state_cache(filename)

    xi = data[:, 0].reshape(ny, nx)
    eta = data[:, 1].reshape(ny, nx)
    x = data[:, 2].reshape(ny, nx)
    y = data[:, 3].reshape(ny, nx)
    w = data[:, 4].reshape(ny, nx)
    psi = data[:, 5].reshape(ny, nx)
    ux = data[:, 6].reshape(ny, nx)
    uy = data[:, 7].reshape(ny, nx)

    prop_map = {
        "vorticity": w,
        "psi": psi,
        "u": ux,
        "v": uy,
    }

    if prop_name not in prop_map:
        raise ValueError(f"Unknown property: {prop_name}")

    prop_data = prop_map[prop_name]
    prop_label = display_prop_name(prop_name)
    meta_label = metadata_text(reynolds_number, nx_meta, ny_meta, u0)

    fig, axes = plt.subplots(1, 2, figsize=figsize)

    ax = axes[0]
    if plot_areas:
        cs = ax.contourf(x, y, prop_data, levels=num_levels, cmap="RdBu_r")
        plt.colorbar(cs, ax=ax, label=prop_label)
    if plot_lines:
        ax.contour(x, y, prop_data, levels=contour_lines, colors="k", alpha=alpha, linewidths=linewidth)
    ax.set_title(f"{prop_label} (x-y domain)")
    ax.set_xlabel("x")
    ax.set_ylabel("y")

    ax = axes[1]
    if plot_areas:
        cs = ax.contourf(xi, eta, prop_data, levels=num_levels, cmap="RdBu_r")
        plt.colorbar(cs, ax=ax, label=prop_label)
    if plot_lines:
        ax.contour(xi, eta, prop_data, levels=contour_lines, colors="k", alpha=alpha, linewidths=linewidth)
    ax.set_title(f"{prop_label} (ξ-η domain)")
    ax.set_xlabel("ξ")
    ax.set_ylabel("η")
    ax.set_aspect("equal")

    if meta_label:
        fig.suptitle(meta_label)

    if meta_label:
        plt.tight_layout(rect=(0.0, 0.0, 1.0, 0.95))
    else:
        plt.tight_layout()

    if save_images:
        try:
            output_path = save_path(prop_name, reynolds_number, nx_meta, ny_meta, u0)
            fig.savefig(output_path, dpi=150, bbox_inches="tight")
            print(f"Saved plot to {output_path}")
        except Exception as e:
            print(f"Failed to save image: {e}")

    if matplotlib.get_backend().lower() != "agg":
        plt.show()
    plt.close(fig)


def main(argv):
    parser = argparse.ArgumentParser(description="Plot state cache with property contours")
    parser.add_argument("filename", help="State cache file")
    parser.add_argument("property", help="Property to plot (vorticity, psi, u, v)")
    parser.add_argument("--levels", type=int, default=NUM_LEVELS, help="Number of contour levels")
    parser.add_argument("--lines", type=int, default=CONTOUR_LINES, help="Number of contour lines")
    parser.add_argument("--width", type=float, default=FIGSIZE_WIDTH, help="Figure width")
    parser.add_argument("--areas", type=int, default=1, help="Plot contour areas (0 or 1)")
    parser.add_argument("--plot-lines", type=int, default=1, help="Plot contour lines (0 or 1)")
    parser.add_argument("--linewidth", type=float, default=0.5, help="Line width for contours")
    parser.add_argument("--alpha", type=float, default=0.3, help="Alpha for contour lines")
    parser.add_argument("--save-images", type=int, default=0, help="Save images (0 or 1)")
    parser.add_argument("--re", type=float, default=None, help="Reynolds number")
    parser.add_argument("--nx", type=int, default=None, help="Grid size in x")
    parser.add_argument("--ny", type=int, default=None, help="Grid size in y")
    parser.add_argument("--u0", type=float, default=None, help="Lid speed u0")

    args = parser.parse_args(argv[1:])
    plot_areas = args.areas != 0
    plot_lines = args.plot_lines != 0
    save_imgs = args.save_images != 0

    plot_property(
        args.filename,
        args.property,
        num_levels=args.levels,
        contour_lines=args.lines,
        figsize_width=args.width,
        plot_areas=plot_areas,
        plot_lines=plot_lines,
        linewidth=args.linewidth,
        alpha=args.alpha,
        save_images=save_imgs,
        reynolds_number=args.re,
        nx_meta=args.nx,
        ny_meta=args.ny,
        u0=args.u0,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))

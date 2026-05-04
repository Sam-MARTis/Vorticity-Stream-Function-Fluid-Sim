import sys
import argparse
from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
import numpy as np

VELOCITY_SCALE = 0.5


def read_centerline_file(path):
	u_lines = []
	v_lines = []
	mode = None
	with open(path, "r", encoding="utf-8") as f:
		for raw_line in f:
			line = raw_line.strip()
			if not line:
				continue
			if line.startswith("#"):
				if line.startswith("# U_CENTER"):
					mode = "u"
				elif line.startswith("# V_CENTER"):
					mode = "v"
				continue
			parts = line.split()
			if len(parts) < 4 or mode not in {"u", "v"}:
				continue
			xi, eta, u_xi, u_eta = map(float, parts[:4])
			if mode == "u":
				u_lines.append((xi, eta, u_xi, u_eta))
			else:
				v_lines.append((xi, eta, u_xi, u_eta))
	return np.asarray(u_lines, dtype=float), np.asarray(v_lines, dtype=float)


def metadata_text(reynolds_number, nx, ny, u0):
	parts = []
	if reynolds_number is not None:
		parts.append(f"Re={reynolds_number:g}")
	if nx is not None and ny is not None:
		parts.append(f"Grid={nx}x{ny}")
	if u0 is not None:
		parts.append(f"u0={u0:g}")
	return " | ".join(parts)


def save_path(reynolds_number, nx, ny, u0):
	images_dir = Path("Images")
	images_dir.mkdir(exist_ok=True)

	suffix_parts = []
	if reynolds_number is not None:
		suffix_parts.append(f"Re{reynolds_number:g}")
	if nx is not None and ny is not None:
		suffix_parts.append(f"N{nx}x{ny}")
	if u0 is not None:
		suffix_parts.append(f"u0{u0:g}")

	suffix = "_" + "_".join(suffix_parts) if suffix_parts else ""
	return images_dir / f"velocity_centerlines{suffix}.png"


def plot_from_file(
	filename,
	scale=1.0,
	linewidth=0.5,
	alpha=0.3,
	save_images=False,
	reynolds_number=None,
	nx=None,
	ny=None,
	u0=None,
):
	u_arr, v_arr = read_centerline_file(filename)
	if u_arr.size == 0 and v_arr.size == 0:
		raise RuntimeError(f"No data read from {filename}")

	fig, ax = plt.subplots(figsize=(6, 6))

	if u_arr.size:
		xi_u = u_arr[:, 0]
		eta_u = u_arr[:, 1]
		u_xi = u_arr[:, 2]
		xi_center = float(xi_u[0])
		xi_plot = xi_center + u_xi * scale
		ax.axvline(x=xi_center, color="k", linestyle="--", linewidth=linewidth, label="_nolegend_")
		ax.plot(xi_plot, eta_u, "-o", color="C0", label="u", linewidth=linewidth)

	if v_arr.size:
		xi_v = v_arr[:, 0]
		eta_v = v_arr[:, 1]
		u_eta = v_arr[:, 3]
		eta_center = float(eta_v[0])
		eta_plot = eta_center + u_eta * scale
		ax.axhline(y=eta_center, color="k", linestyle="--", linewidth=linewidth, label="_nolegend_")
		ax.plot(xi_v, eta_plot, "-s", color="C1", label="v", linewidth=linewidth)

	ax.set_xlabel("ξ")
	ax.set_ylabel("η")
	ax.set_aspect("equal", adjustable="box")
	ax.grid(True, linestyle=":", alpha=alpha)
	ax.legend()

	ax.relim()
	ax.autoscale_view()
	xmin, xmax = ax.get_xlim()
	ymin, ymax = ax.get_ylim()
	half = 0.55 * max(xmax - xmin, ymax - ymin)
	cx = 0.5 * (xmin + xmax)
	cy = 0.5 * (ymin + ymax)
	ax.set_xlim(cx - half, cx + half)
	ax.set_ylim(cy - half, cy + half)

	# Add secondary axes showing velocities scaled from coordinate displacements
	# Top axis: horizontal displacement -> u velocity (if u data present)
	if u_arr.size:
		def _xi_to_u(x):
			return (x - xi_center) / scale

		def _u_to_x(u):
			return xi_center + u * scale

		sec_x = ax.secondary_xaxis("top", functions=(_xi_to_u, _u_to_x))
		sec_x.set_xlabel("u (velocity)")

	# Right axis: vertical displacement -> v velocity (if v data present)
	if v_arr.size:
		def _eta_to_v(y):
			return (y - eta_center) / scale

		def _v_to_eta(vv):
			return eta_center + vv * scale

		sec_y = ax.secondary_yaxis("right", functions=(_eta_to_v, _v_to_eta))
		sec_y.set_ylabel("v (velocity)")

	meta_label = metadata_text(reynolds_number, nx, ny, u0)
	if meta_label:
		fig.suptitle(meta_label)
		plt.tight_layout(rect=(0.0, 0.0, 1.0, 0.95))
	else:
		plt.tight_layout()

	if save_images:
		try:
			output_path = save_path(reynolds_number, nx, ny, u0)
			fig.savefig(output_path, dpi=150, bbox_inches="tight")
			print(f"Saved plot to {output_path}")
		except Exception as e:
			print(f"Failed to save image: {e}")

	if matplotlib.get_backend().lower() != "agg":
		plt.show()
	plt.close(fig)


def main(argv):
	parser = argparse.ArgumentParser(description="Plot centerline velocities")
	parser.add_argument("filename", help="Centerline data file")
	parser.add_argument("--linewidth", type=float, default=0.5, help="Line width for plots")
	parser.add_argument("--alpha", type=float, default=0.3, help="Alpha for grid")
	parser.add_argument("--save-images", type=int, default=0, help="Save images (0 or 1)")
	parser.add_argument("--re", type=float, default=None, help="Reynolds number")
	parser.add_argument("--nx", type=int, default=None, help="Grid size in x")
	parser.add_argument("--ny", type=int, default=None, help="Grid size in y")
	parser.add_argument("--u0", type=float, default=None, help="Lid speed u0")
	
	args = parser.parse_args(argv[1:])
	save_imgs = args.save_images != 0
	plot_from_file(
		args.filename,
		scale=VELOCITY_SCALE,
		linewidth=args.linewidth,
		alpha=args.alpha,
		save_images=save_imgs,
		reynolds_number=args.re,
		nx=args.nx,
		ny=args.ny,
		u0=args.u0,
	)
	return 0


if __name__ == "__main__":
	raise SystemExit(main(sys.argv))
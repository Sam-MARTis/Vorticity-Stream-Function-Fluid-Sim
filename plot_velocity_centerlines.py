import sys

import matplotlib
import matplotlib.pyplot as plt
import numpy as np

VELOCITY_SCALE = 1.0


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


def plot_from_file(filename, scale=1.0):
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
		ax.axvline(x=xi_center, color="k", linestyle="--", linewidth=0.8, label="_nolegend_")
		ax.plot(xi_plot, eta_u, "-o", color="C0", label="u")

	if v_arr.size:
		xi_v = v_arr[:, 0]
		eta_v = v_arr[:, 1]
		u_eta = v_arr[:, 3]
		eta_center = float(eta_v[0])
		eta_plot = eta_center + u_eta * scale
		ax.axhline(y=eta_center, color="k", linestyle="--", linewidth=0.8, label="_nolegend_")
		ax.plot(xi_v, eta_plot, "-s", color="C1", label="v")

	ax.set_xlabel("xi")
	ax.set_ylabel("eta")
	ax.set_aspect("equal", adjustable="box")
	ax.grid(True, linestyle=":", alpha=0.4)
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

	if matplotlib.get_backend().lower() != "agg":
		plt.show()
	plt.close(fig)


def main(argv):
	if len(argv) < 2:
		return 1
	plot_from_file(argv[1], scale=VELOCITY_SCALE)
	return 0


if __name__ == "__main__":
	raise SystemExit(main(sys.argv))
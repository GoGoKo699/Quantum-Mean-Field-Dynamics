"""Render the recorded coupling experiment; no new simulation is performed."""
from pathlib import Path
import json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
data = json.loads((ROOT / "verification/random_time_coupling_check.json").read_text())
plt.rcParams.update({"svg.fonttype": "none", "svg.hashsalt": "mean-field-coupling",
                     "font.size": 11, "axes.spines.top": False,
                     "axes.spines.right": False})
fig, ax = plt.subplots(figsize=(7.4, 4.6), layout="constrained")
for width in sorted({r["H"] for r in data["records"]}, reverse=True):
    rows = sorted((r for r in data["records"] if r["H"] == width),
                  key=lambda r: r["h"])
    ax.loglog([r["h"] for r in rows], [r["std"] / width for r in rows],
              marker="o", linewidth=1.8, markersize=4.5, label=f"H = {width:g}")
h = np.array([min(r["h"] for r in data["records"]),
              max(r["h"] for r in data["records"])])
ax.loglog(h, 0.1*h, "--", color="0.3", linewidth=1.4,
          label="Slope 1 reference")
ax.set_xlabel("Euler mesh h")
ax.set_ylabel("Standard deviation of correction / H")
ax.set_title("Random-time Dynkin correction scaling", loc="left", pad=12)
ax.grid(which="both", alpha=0.18)
ax.legend(frameon=False, fontsize=9, ncol=2)
destination = ROOT / "assets/coupling-scaling.svg"
destination.parent.mkdir(parents=True, exist_ok=True)
fig.savefig(destination, metadata={"Date": None})
print(destination)

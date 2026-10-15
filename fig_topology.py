"""fig_topology.py - draws the 5-router worked example used in Section 3."""
import json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Circle

d = json.load(open("worked.json"))
pos = {"A": (0.0, 0.0), "B": (1.2, 0.55), "C": (2.4, 0.55), "D": (3.6, 0.55), "E": (4.8, 0.0)}
rad = {("A", "E"): 0.32, ("B", "D"): -0.45}
lab = {("A", "E"): (2.4, -0.62), ("B", "D"): (2.4, 1.38)}

fig, ax = plt.subplots(figsize=(6.4, 2.6))
ax.set_xlim(-0.5, 5.3); ax.set_ylim(-0.95, 1.7); ax.axis("off"); ax.set_aspect("equal")
for a, b, w in d["edges"]:
    r = rad.get((a, b), 0.0)
    ax.add_patch(FancyArrowPatch(pos[a], pos[b], arrowstyle="-", lw=1.8,
                                 color="#44546a", connectionstyle=f"arc3,rad={r}",
                                 shrinkA=14, shrinkB=14))
    if (a, b) in lab:
        x, y = lab[(a, b)]
    else:
        x = (pos[a][0] + pos[b][0]) / 2; y = (pos[a][1] + pos[b][1]) / 2 + 0.17
    ax.text(x, y, str(w), ha="center", va="center", fontsize=10, color="#b03a2e", weight="bold",
            bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none"))
for n, (x, y) in pos.items():
    ax.add_patch(Circle((x, y), 0.2, fc="#e8f0f7", ec="#1f3b57", lw=1.8, zorder=3))
    ax.text(x, y, n, ha="center", va="center", fontsize=11, weight="bold", zorder=4)
fig.tight_layout(pad=0.3)
fig.savefig("fig_topology.png", dpi=170)
print("ok")

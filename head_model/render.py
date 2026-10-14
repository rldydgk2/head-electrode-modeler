"""Deterministic surface previews; no geometry generation or electric-field solve."""

from pathlib import Path
import numpy as np

from .poly import read_poly


def render(path, output, *, title="Head model", max_faces=300000):
    """Render two views of the actual skin and electrode surface facets."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection

    surface = read_poly(path)
    skin = np.flatnonzero(surface.markers == 1005)
    other = np.flatnonzero((surface.markers < 1000))
    if not len(skin):
        skin = np.arange(len(surface.triangles))
        other = np.array([], dtype=int)
    skin = skin[:: max(1, int(np.ceil(len(skin) / max_faces)))]
    fig = plt.figure(figsize=(14, 8), facecolor="#101b2a")
    fig.suptitle(title, color="white", fontsize=25, x=0.06, ha="left", y=0.94)
    fig.text(
        0.06,
        0.875,
        "TETGEN HEAD MODEL  /  SURFACE GEOMETRY  /  UNITS: mm",
        color="#91adbf",
        fontsize=10,
    )
    for panel, azimuth in enumerate((-55, 125), 1):
        ax = fig.add_subplot(1, 2, panel, projection="3d", facecolor="#101b2a")
        selected = np.concatenate((skin, other))
        colors = ["#a8c8d4"] * len(skin) + [
            "#efb65b" if marker == -1 else "#ed785e"
            for marker in surface.markers[other]
        ]
        ax.add_collection3d(
            Poly3DCollection(
                surface.points[surface.triangles[selected]],
                facecolors=colors,
                shade=True,
                linewidth=0,
                alpha=1,
            )
        )
        bounds = surface.points[surface.triangles[skin]].reshape(-1, 3)
        low, high = bounds.min(axis=0), bounds.max(axis=0)
        center, radius = (low + high) / 2, max(high - low) / 2 + 12
        ax.set(
            xlim=(center[0] - radius, center[0] + radius),
            ylim=(center[1] - radius, center[1] + radius),
            zlim=(center[2] - radius, center[2] + radius),
        )
        ax.set_box_aspect((1, 1, 1))
        ax.view_init(elev=15, azim=azimuth)
        ax.set_axis_off()
    fig.text(
        0.06,
        0.07,
        "Skin surface" + ("     /     Electrode assembly" if len(other) else ""),
        color="#d2e0e9",
        fontsize=12,
    )
    fig.text(
        0.06,
        0.035,
        "Surface preview only. Internal tissues and tumor are hidden; no field solution is shown.",
        color="#91adbf",
        fontsize=9,
    )
    fig.subplots_adjust(left=0, right=1, bottom=0.09, top=0.88, wspace=-0.12)
    Path(output).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=150, facecolor=fig.get_facecolor())
    plt.close(fig)

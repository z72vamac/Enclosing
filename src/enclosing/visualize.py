"""Visualizacion matplotlib."""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Rectangle, Polygon as MplPolygon


def plot_solution(poly, centers, R, demands=None, save=None, show=False,
                  title=None, shape="circle"):
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.set_aspect("equal")
    # poligono
    x, y = poly.exterior.xy
    ax.fill(x, y, alpha=0.15, label="region")
    ax.plot(x, y, "k-", lw=1.2)
    for ring in poly.interiors:
        xi, yi = ring.xy
        ax.fill(xi, yi, color="white")
        ax.plot(xi, yi, "k-", lw=1.2)
    # piezas
    for (cx, cy) in np.asarray(centers):
        if shape == "square":
            c = Rectangle((cx - R, cy - R), 2 * R, 2 * R, alpha=0.25,
                          edgecolor="C0", facecolor="C0")
        elif shape == "circle":
            c = Circle((cx, cy), R, alpha=0.25, edgecolor="C0", facecolor="C0")
        else:
            raise ValueError(f"shape desconocido: {shape!r}")
        ax.add_patch(c)
        ax.plot(cx, cy, "C0o", ms=4)
    if demands is not None and len(demands):
        ax.plot(demands[:, 0], demands[:, 1], "k.", ms=2, alpha=0.4, label="demanda")
    ax.autoscale_view()
    if title:
        ax.set_title(title)
    ax.legend(loc="best")
    plt.tight_layout()
    if save:
        fig.savefig(save, dpi=150)
    if show:
        plt.show()
    plt.close(fig)
    return save


def plot_tour(poly, centers, tour, R=None, save=None, show=False, title=None,
              shape="circle"):
    C = np.asarray(centers, dtype=float)
    t = list(tour) + [tour[0]]
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.set_aspect("equal")
    x, y = poly.exterior.xy
    ax.fill(x, y, alpha=0.12)
    ax.plot(x, y, "k-", lw=1.2)
    if R is not None:
        for (cx, cy) in C:
            if shape == "square":
                ax.add_patch(Rectangle((cx - R, cy - R), 2 * R, 2 * R,
                                       alpha=0.12, edgecolor="C0",
                                       facecolor="C0"))
            else:
                ax.add_patch(Circle((cx, cy), R, alpha=0.12,
                                    edgecolor="C0", facecolor="C0"))
    ax.plot(C[t, 0], C[t, 1], "C1-", lw=1.0, alpha=0.9, label="tour")
    ax.plot(C[:, 0], C[:, 1], "C0o", ms=3)
    ax.plot(C[t[0], 0], C[t[0], 1], "rs", ms=5, label="inicio")
    if title:
        ax.set_title(title)
    ax.legend(loc="best")
    plt.tight_layout()
    if save:
        fig.savefig(save, dpi=150)
    if show:
        plt.show()
    plt.close(fig)
    return save

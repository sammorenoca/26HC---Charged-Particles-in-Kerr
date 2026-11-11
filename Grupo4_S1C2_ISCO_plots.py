import csv
from collections import OrderedDict

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import AutoMinorLocator

plt.rcParams["font.family"] = "serif"
plt.rcParams["font.serif"] = ["Times New Roman", "Liberation Serif", "STIXGeneral", "DejaVu Serif"]
plt.rcParams["mathtext.fontset"] = "stix"

LAB = 17   
TICK = 15  


def read(archivo="curvas.csv"):
    datos = OrderedDict()
    with open(archivo) as f:
        for fila in csv.DictReader(f):
            clave = (fila["panel"], float(fila["valor"]))
            datos.setdefault(clave, []).append((float(fila["r"]), float(fila["V"])))
    out = OrderedDict()
    for (panel, valor), pts in datos.items():
        out.setdefault(panel, []).append((valor, np.array(pts)))
    return out


def style(ax, xlabel=r"$r$", ylabel=r"$V_{\rm eff}$"):
    ax.set_box_aspect(1)
    ax.set_xlabel(xlabel, fontsize=LAB)
    ax.set_ylabel(ylabel, fontsize=LAB)
    ax.tick_params(which="both", direction="in", top=True, right=True, labelsize=TICK)
    ax.tick_params(which="major", length=6, width=1.0)
    ax.tick_params(which="minor", length=3, width=1.0)
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.grid(which="major", linestyle="--", alpha=0.45)
    ax.grid(which="minor", linestyle=":", alpha=0.25)
    for s in ax.spines.values():
        s.set_visible(True)
        s.set_color("black")
        s.set_linewidth(1.0)


def plot(ax, curvas, etiqueta):
    n = len(curvas)
    colores = cm.viridis(np.linspace(0.1, 0.9, n))
    for (valor, pts), c in zip(curvas, colores):
        ax.plot(pts[:, 0], pts[:, 1], "-", linewidth=2.0, color=c, label=etiqueta(valor))
    ax.legend(frameon=False, fontsize=13)


def exp_tex(x, base=1e-3, nombre=r"B^*"):
    return rf"${nombre}={x/base:g}\times10^{{{int(np.log10(base))}}}$"


def guardar(fig, nombre):
    fig.set_layout_engine("constrained")
    fig.savefig(nombre, dpi=600, bbox_inches="tight", pad_inches=0.1)
    plt.close(fig)


datos = read()

fig, ax = plt.subplots(figsize=(6, 6))
plot(ax, datos["fig1"], lambda v: exp_tex(v, 1e-3, r"B^*"))
style(ax)
ax.set_xlim(0, 20)
ax.set_ylim(0.930, 0.960)
guardar(fig, "fig1_potencial_B.pdf")

fig, ax = plt.subplots(figsize=(6, 6))
plot(ax, datos["fig2"], lambda v: rf"$Q^*={v/1e-4:g}\times10^{{-4}}$")
style(ax)
ax.set_xlim(3, 10)
ax.set_ylim(0.943, 0.947)
guardar(fig, "fig2_potencial_Q.pdf")

fig, ax = plt.subplots(figsize=(6, 6))
plot(ax, datos["fig3"], lambda v: rf"$a={v:g}$")
style(ax)
ax.set_xlim(0, 15)
ax.set_ylim(0.936, 0.950)
guardar(fig, "fig3_potencial_a.pdf")

fig, axs = plt.subplots(1, 3, figsize=(18, 6))
paneles = [("fig4a", 0.0, (0, 50), (0.80, 1.15)),
           ("fig4b", 0.4, (0, 50), (0.80, 1.15)),
           ("fig4c", 0.6, (0, 35), (0.85, 1.00))]

for ax, (clave, a, xl, yl), letra in zip(axs, paneles, "abc"):
    plot(ax, datos[clave], lambda v: rf"$L={f'{v:.4f}'.rstrip('0').rstrip('.')}$")
    style(ax)
    ax.set_xlim(*xl)
    ax.set_ylim(*yl)
    ax.text(-0.02, 1.04, f"({letra})", transform=ax.transAxes, fontsize=LAB, ha="right", va="bottom")
    ax.text(0.97, 0.05, rf"$a={a:g}$", transform=ax.transAxes, fontsize=14, ha="right", va="bottom")

guardar(fig, "fig4_potencial_L.pdf")
print("Figuras escritas: fig1_potencial_B.pdf, fig2_potencial_Q.pdf, fig3_potencial_a.pdf, fig4_potencial_L.pdf")
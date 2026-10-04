"""Recompute the circular-ring enrichment factors with the correct denominator.

The original analysis divided the local share of the *neighbour* class by the
global share of the *changed-to* class. The enrichment factor (Verburg et al.,
2004) divides by the global share of the *neighbour* class:

    F[i, k, d] = (n_k,d / n_d) / (N_k / N)

Within one panel the wrong denominator is a single constant, so the correction
is exact and needs no raster data:

    log F_correct = log F_original + ln(N_changed_to / N_neighbour)

All inputs (per-ring log F values and per-class pixel counts) are read from the
outputs saved in the original notebook.
"""
import ast
import json
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).parent
NOTEBOOK = HERE.parent / "original" / "land_use_enrichment_factor_2024.ipynb"
CLASSES = {"u": "Urban", "i": "Industrial", "r": "Rural", "a": "Agricultural"}


def saved_output(nb, source):
    """Return the text output of the code cell whose source is exactly `source`."""
    for cell in nb["cells"]:
        if cell["cell_type"] == "code" and "".join(cell["source"]).strip() == source:
            for out in cell["outputs"]:
                text = out.get("text") or out.get("data", {}).get("text/plain", "")
                if text:
                    return "".join(text)
    raise LookupError(source)


def load_inputs():
    nb = json.loads(NOTEBOOK.read_text(encoding="utf-8"))

    # per-class pixel counts in the 1390 map ("Image 2" table)
    table = saved_output(nb, nb_source_counts(nb)).split("Image 2:")[1]
    counts = {name: int(re.search(rf"{name}\s+(\d+)", table).group(1)) for name in CLASSES.values()}

    # per-ring log enrichment factors (circular rings)
    raw = saved_output(nb, "log_enrich_f_circ_dict.items()")
    raw = raw.replace("dict_items(", "", 1).rstrip()[:-1].replace("-inf", "None")
    log_f = {k: {r: (-np.inf if v is None else v) for r, v in d.items()}
             for k, d in ast.literal_eval(raw)}
    return counts, log_f


def nb_source_counts(nb):
    for cell in nb["cells"]:
        src = "".join(cell["source"]).strip()
        if cell["cell_type"] == "code" and 'print("\\nImage 2:")' in src:
            return src
    raise LookupError("pixel-count cell")


def correct(counts, log_f):
    corrected = {}
    for key, rings in log_f.items():
        changed_to, neighbour = CLASSES[key[0]], CLASSES[key[1]]
        shift = np.log(counts[changed_to] / counts[neighbour])
        corrected[key] = {r: v + shift for r, v in rings.items()}
    return corrected


def plot(log_f, corrected, path):
    fig, axs = plt.subplots(4, 4, figsize=(16, 15), sharex=True, sharey=True)
    for ax, key in zip(axs.flat, log_f):
        rings = sorted(log_f[key])
        orig = np.array([log_f[key][r] for r in rings])
        corr = np.array([corrected[key][r] for r in rings])
        ax.axhline(0, color="black", linewidth=0.8)
        ax.plot(rings, orig, "o--", color="gray", label="Original (2024)")
        ax.plot(rings, corr, "o-", color="tab:blue", label="Corrected denominator")
        a, b = CLASSES[key[0]], CLASSES[key[1]]
        ax.set_title(f"{key}: {b} around new {a}", fontsize=10)
        ax.grid(True, alpha=0.4)
    for ax in axs[-1]:
        ax.set_xlabel("Ring")
    for ax in axs[:, 0]:
        ax.set_ylabel("log(F)")
    axs[0, 0].legend(fontsize=9, loc="lower left")
    fig.suptitle("Log enrichment factor, circular rings: original vs corrected denominator\n"
                 "(above 0 = attraction, below 0 = repulsion; diagonal panels are unchanged)")
    plt.tight_layout()
    plt.savefig(path, dpi=110)


if __name__ == "__main__":
    counts, log_f = load_inputs()
    corrected = correct(counts, log_f)
    print("Pixel counts (1390):", counts)
    print(f"{'pair':6s} {'ring 1 orig':>12s} {'ring 1 corr':>12s} {'ring 8 orig':>12s} {'ring 8 corr':>12s}")
    for k in log_f:
        print(f"{k:6s} {log_f[k][1]:12.2f} {corrected[k][1]:12.2f} {log_f[k][8]:12.2f} {corrected[k][8]:12.2f}")
    plot(log_f, corrected, HERE.parent / "figures" / "01_corrected_vs_original.png")

"""Part 2: recompute the network-distance enrichment factors of the 50 sample parcels.

The 2024 notebook (original/part2_network/) measured, for each sample parcel, the
area of every land use reachable along the street network within 0-50, 50-100,
100-150 and 150-200 m, and pooled these areas per land use of the central parcel.

Two things are corrected here, exactly, from the numbers saved in that notebook:

1. Denominator. The original divided by the land-use composition of the pooled
   neighbourhoods themselves (and mixed in vacant land in the total). The
   enrichment factor compares against the whole study area:
       F = (A_k,band / A_band) / (A_k,study / A_study)      (vacant 'bayer' excluded)
2. Self-inclusion. Every sample parcel was counted inside its own 0-50 m band,
   which inflates same-class values. Its own area is removed from that band.

Only the 50 sample parcels can be recomputed; the original rasters/vectors are
not available. Results are therefore a small-sample illustration.
"""
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).parent
NOTEBOOK = HERE.parent / "original" / "part2_network" / "network_enrichment_factor_2024.ipynb"

BANDS = ["0-50 meters", "50-100 meters", "100-150 meters", "150-200 meters"]
BAND_LABELS = ["0–50", "50–100", "100–150", "150–200"]
CLASSES = ["maskooni", "tejari", "khadamati", "park", "sanati"]
NAMES = {"maskooni": "Residential", "tejari": "Commercial", "khadamati": "Services",
         "park": "Park", "sanati": "Industrial"}
CENTRAL = ["maskooni", "tejari", "park", "khadamati"]   # sample classes with saved results


def outputs(nb):
    """Map each code cell's printed text to its source, for lookup by content."""
    for cell in nb["cells"]:
        if cell["cell_type"] != "code":
            continue
        for out in cell.get("outputs", []):
            text = out.get("text") or out.get("data", {}).get("text/plain", "")
            if text:
                yield "".join(cell["source"]), "".join(text)


def first_json_after(text, marker):
    start = text.index(marker) + len(marker)
    return json.JSONDecoder().raw_decode(text[start:].lstrip())[0]


def load(nb):
    band_areas, original_ef, samples, study = {}, {}, None, None
    for src, text in outputs(nb):
        for c in CENTRAL:
            prefix = "" if c == "maskooni" else f"{c}_"
            if text.startswith("Area of each label") and f"{prefix}area_per_label_per_range" in src \
                    and (c != "maskooni" or "tejari_" not in src and "park_" not in src and "khadamati_" not in src):
                band_areas.setdefault(c, first_json_after(text, "range:"))
            if text.startswith("Enrichment factors:") and f"{prefix}enrichment_factors[range_key][label] =" in src:
                original_ef.setdefault(c, first_json_after(text, "Enrichment factors:"))  # first match
        if "json.dumps(sample_50_centroids_data" in src and samples is None:
            samples = json.loads(text)
        if text.startswith("Total area of the study area"):
            study = eval(text.split("study area: ", 2)[2].strip())   # printed dict literal
    return band_areas, original_ef, samples, study


def correct(band_areas, samples, study):
    study_total = sum(v for k, v in study.items() if k != "bayer")
    share = {k: study[k] / study_total for k in CLASSES}
    own_area = {c: sum(s["attributes"]["SHAPE_Area"] for s in samples.values()
                       if s["attributes"]["landuse2"] == c) for c in CENTRAL}
    n = {c: sum(s["attributes"]["landuse2"] == c for s in samples.values()) for c in CENTRAL}

    corrected = {}
    for c in CENTRAL:
        areas = {b: dict(band_areas[c][b]) for b in BANDS}
        areas["0-50 meters"][c] -= own_area[c]          # remove the parcel itself
        corrected[c] = {}
        for b in BANDS:
            band_total = sum(areas[b][k] for k in CLASSES)
            corrected[c][b] = {k: (areas[b][k] / band_total) / share[k] if band_total else np.nan
                               for k in CLASSES}
    return corrected, n, share


def plot(original_ef, corrected, n, path):
    fig, axs = plt.subplots(len(CENTRAL), len(CLASSES), figsize=(18, 13), sharex=True)
    x = np.arange(len(BANDS))
    for i, c in enumerate(CENTRAL):
        for j, k in enumerate(CLASSES):
            ax = axs[i, j]
            orig = [original_ef[c][b][k] for b in BANDS]
            corr = [corrected[c][b][k] for b in BANDS]
            ax.axhline(1, color="black", linewidth=0.8)
            ax.plot(x, orig, "o--", color="gray", label="Original (2024)")
            ax.plot(x, corr, "o-", color="tab:blue", label="Corrected")
            ax.set_title(f"{NAMES[k]} around {NAMES[c]} (n={n[c]})", fontsize=10)
            ax.set_xticks(x, BAND_LABELS)
            ax.grid(True, alpha=0.4)
            if j == 0:
                ax.set_ylabel("Enrichment factor F")
    for ax in axs[-1]:
        ax.set_xlabel("Network distance band (m)")
    axs[0, 0].legend(fontsize=9)
    fig.suptitle("Network-distance enrichment factors for the 50 sample parcels: original vs corrected\n"
                 "(F > 1 = over-represented, F < 1 = under-represented; n = sample parcels of that class)")
    plt.tight_layout()
    plt.savefig(path, dpi=110)


if __name__ == "__main__":
    nb = json.loads(NOTEBOOK.read_text(encoding="utf-8"))
    band_areas, original_ef, samples, study = load(nb)
    assert set(band_areas) == set(original_ef) == set(CENTRAL), (band_areas.keys(), original_ef.keys())
    corrected, n, share = correct(band_areas, samples, study)
    print("Sample parcels per class:", n)
    print("Study-area shares (excl. vacant):", {k: round(v, 4) for k, v in share.items()})
    for c in CENTRAL:
        print(f"\n{NAMES[c]} (n={n[c]})   band: " + "  ".join(f"{k[:6]:>13s}" for k in CLASSES))
        for b in BANDS:
            print(f"  {b:15s}" + "  ".join(f"{original_ef[c][b][k]:5.2f}->{corrected[c][b][k]:5.2f}" for k in CLASSES))
    plot(original_ef, corrected, n, HERE.parent / "figures" / "02_part2_network_ef_corrected.png")

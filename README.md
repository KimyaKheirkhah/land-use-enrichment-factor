# Land-Use Neighbourhood Analysis with the Enrichment Factor

Measuring how land-use change is shaped by its surroundings. For every pixel that changed land use between two maps, this analysis asks which land uses were over- or under-represented around it, ring by ring, out to eight rings. The resulting attraction/repulsion curves are the empirical basis for the neighbourhood rules of cellular-automaton land-use models.

> Developed in **June 2024** for the *Spatial & Land Use Planning* course of my MSc in GIS (K. N. Toosi University of Technology).
> The notebook is published as originally written, with its explanations in Persian. This README describes it in English and adds a [2026 review](#review-2026) with one exact correction.
>
> **Part 1** (below) works on raster land-use maps with straight-line rings. **[Part 2](#part-2-network-distance-neighbourhoods-on-parcels)** takes the same idea to vector parcels, measuring distance along the street network.

---

## The problem

Land-use change is not random in space. New urban areas appear next to existing ones, industry clusters, and farmland gives way at the urban fringe. To model this, for example with a cellular automaton (see my [urban growth CA project](https://github.com/KimyaKheirkhah/urban-growth-cellular-automata)), we need to know **how strongly each land use attracts or repels each other land use, and how that effect fades with distance**.

The **enrichment factor** (Verburg et al., 2004) quantifies this:

$$F_{k,d} = \frac{n_{k,d} / n_d}{N_k / N}$$

where $n_{k,d}/n_d$ is the share of land use *k* in the ring at distance *d* around a location, and $N_k/N$ is the share of *k* in the whole map.

- **F > 1** (log F > 0): *k* is over-represented nearby, i.e. **attraction**.
- **F < 1** (log F < 0): *k* is under-represented, i.e. **repulsion**.

## Data and approach

Two classified land-use maps of the same area, for **1380 and 1390** in the Iranian calendar (about 2001 and 2011), with urban, industrial, rural, agricultural and barren classes.

![Land-use maps](figures/orig_01_land_use_maps.png)

1. **Change detection.** For each target class (urban, industrial, rural, agricultural), find the pixels that changed *into* that class between the two dates.

   ![Changed pixels](figures/orig_02_changed_pixels.png)

2. **Neighbourhood rings.** Around every changed pixel, count each land use in 8 concentric rings. Two ring shapes were implemented and compared: **square rings** (3×3, 5×5, … 17×17 windows minus the inner window) and **circular rings** (pixels at distance *r ± 0.5*).

   ![Circular rings](figures/orig_04_circular_rings.png)

3. **Enrichment factors.** Compute F for every changed pixel and ring, then average per ring, for all **16 combinations** of changed-to class × neighbour class (urban–urban, urban–industrial, …, agricultural–agricultural).

   *Square rings:*
   ![EF square rings](figures/orig_03_ef_square_rings.png)

   *Circular rings, log scale:*
   ![Log EF circular rings](figures/orig_05_log_ef_circular_rings.png)

4. **Completing the curves.** Where a ring contains none of a land use, F = 0 and log F = −∞, so the curve starts only at a later ring. As specified in the course, these near rings were filled by fitting a logarithmic function (*a·ln(d) + b*) to the available rings and extrapolating it back to the first ring, so that every curve reaches the axis.

   ![Fitted curves](figures/orig_06_log_ef_fitted.png)

**What the same-class curves show** (urban–urban, industrial–industrial, rural–rural, agricultural–agricultural): new development is strongly attracted to existing land of the same type. The effect is strongest in the first ring and decays steadily with distance, which is the classic pattern behind neighbourhood rules in urban CA models. The strongest clustering is for rural (log F ≈ 5.4 in ring 1) and industrial land (≈ 4.1).

## Review (2026)

### The denominator: an exact correction

In the original code, the global share in the denominator was the share of the **changed-to** class, not of the **neighbour** class being counted. For example, for agricultural land around new urban pixels (`ua`), the agricultural share in the rings was divided by the *urban* share of the map (2.2%) instead of the *agricultural* share (32%).

This leaves the four same-class panels unchanged, but shifts the other 12, in both the square-ring and the circular-ring results. Because the wrong denominator is a single constant per panel, the correction is **exact** and needs no access to the original rasters:

$$\log F_{\text{correct}} = \log F_{\text{original}} + \ln\frac{N_{\text{changed-to}}}{N_{\text{neighbour}}}$$

[`corrected/correct_denominator.py`](corrected/correct_denominator.py) reads the per-ring values and the pixel counts saved in the original notebook and applies this correction. Only the circular-ring values were printed in full, so the corrected figure covers the circular rings:

![Corrected vs original](figures/01_corrected_vs_original.png)

The shape of every curve is unchanged, but several interpretations change:

| Pair | Original reading | Corrected reading |
|---|---|---|
| `ua`: agricultural around new urban | strong attraction (log F ≈ +1.7 → +2.8) | slight repulsion close by, **neutral** further out (−1.0 → +0.2) |
| `ia`: agricultural around new industrial | strong attraction | mild repulsion, approaching neutral |
| `ra`: agricultural around new rural | very strong attraction (≈ +4 → +5) | **close to neutral** (−0.6 → +0.6) |
| `ar`: rural around new agricultural | very strong repulsion (≈ −7 → −4.5) | repulsion close by, **neutral** from ring 5 |
| `ir`: rural around new industrial | repulsion | **attraction** from ring 4 onward |
| `ru`: urban around new rural | attraction at rings 7–8 | **repulsion** at all available rings |

### What I would do differently

- **Measure neighbourhoods on the earlier map.** Neighbour counts were taken from the 1390 map, *after* the change, so they include the converted pixels themselves and other conversions from the same period. The usual approach is the 1380 map: what surrounded a location *before* it changed. This inflates the same-class curves somewhat, and it can't be corrected without the original rasters.
- **Treat F = 0 as information.** A ring with none of a land use is evidence of strong repulsion, not missing data. Extrapolating a curve over it (step 4) produces a smooth, complete weight function for a CA model, but the extrapolated values should be labelled as modelled, not observed. With the data, I would also report how many changed pixels each ring is based on, since some averages rest on very few pixels.
- **Handle map edges.** Rings that extend past the map edge are only partially observed, but are divided by the full ring size. Dividing by the number of *valid* pixels in each ring removes this bias.
- **Vectorize.** The nested pixel loops can be replaced by a convolution per ring and per class, which gives identical counts in seconds.
- **Close the loop.** Use the corrected curves to set the neighbourhood weights of a CA model (as in constrained CA, White & Engelen, 1993) and validate the simulation of 1380 → 1390 against the observed 1390 map.

## Part 2: Network-distance neighbourhoods on parcels

The second assignment of the course moved the same question from pixels to **real parcels**, and from straight-line rings to **distance along the street network**: what a resident or shop can actually reach on foot or by car.

### Data and approach

A vector dataset of an urban district (projected coordinates, metres): land-use parcels classified as residential (*maskooni*), commercial (*tejari*), services (*khadamati*), park, industrial (*sanati*) and vacant (*bayer*), plus restricted zones and a two-level street network.

![Parcels, restricted zones and streets](figures/part2_orig_01_land_use_streets.png)

1. **Street graph.** Every street segment became an edge of a `networkx` graph, weighted by its length. All-pairs Dijkstra gave an origin–destination cost matrix between street nodes.
2. **Parcels to the network.** Each parcel's centroid was linked to its nearest street node with a KD-tree.

   ![Parcels and their nearest street nodes](figures/part2_orig_02_nearest_nodes.png)

3. **Network service areas.** From a parcel's street node, Dijkstra with distance cut-offs gave the nodes reachable within 0–50, 50–100, 100–150 and 150–200 m. A parcel belongs to a band when its own nearest node falls in that band.

   ![Network distance bands around one parcel](figures/part2_orig_03_service_area_bands.png)

4. **Area-weighted enrichment factors.** For a sample of 50 parcels, the area of each land use in each band was pooled per land use of the central parcel and turned into enrichment factors, now weighted by **parcel area** instead of pixel counts.
5. Several further steps were explored: mapping a "neighbourhood effect" back onto every parcel, extending the analysis to all 5,797 parcels, and clustering the land uses by their enrichment profiles (distance matrix, dendrogram, MDS).

### Review (2026)

Steps 1–3 are sound, and network distance is a genuine improvement over straight-line rings for urban land use. Parts of the rest did not work, and the notebook is published as the exploratory working file it was:

- **The final effect maps carry no information.** In the "neighbourhood effect" maps, every parcel receives the same value per panel; a single number per land use was assigned to all parcels instead of computing an effect for each parcel. One attempt shows that all 2,907 values are 0.

  ![Uniform effect maps](figures/part2_orig_04_uniform_effect_maps.png)

- **Several different quantities are called "enrichment factor".** Some later cells compute the share of a land use's area within a buffer around *any street*, which is high for almost every class in a dense street grid and says little about neighbours. Those are not enrichment factors.
- **The denominator came from the neighbourhoods, not the study area.** The "global" share was computed from the pooled neighbourhoods themselves, with vacant land included in the total but excluded from the bands.
- **Each sample parcel was counted in its own 0–50 m band.** This inflates every same-class value at the shortest distance.
- **The sample is the first 50 rows of the table**, clustered in the east of the district: 29 residential, 13 commercial, 5 park, 1 services, 2 vacant and no industrial parcel. The full run over all parcels was started but not completed.

### Recomputed sample enrichment factors

The notebook saved the pooled area of every land use in every band, the area of each sample parcel, and the whole-study-area totals. That is enough to recompute the sample enrichment factors **exactly**, with the study-area denominator and without each parcel counting itself:

$$F_{k,b} = \frac{A_{k,b} / A_b}{A_k / A} \qquad \text{(vacant land excluded throughout)}$$

[`corrected/part2_network_sample_ef.py`](corrected/part2_network_sample_ef.py) reads those saved outputs and draws the comparison:

![Part 2: original vs corrected](figures/02_part2_network_ef_corrected.png)

Read this as a **small-sample illustration**, not a result for the district. Panels based on 5 parks or a single services parcel are anecdotes, and industrial land is only 0.1% of the area, so its values are unstable. With that caveat, the corrections change some readings:

- **Commercial clusters with commercial** (F ≈ 2.9 in the first 50 m), and this survives removing each parcel from its own band, so it is not an artefact.
- **"Parks attract parks" was an artefact.** The original F ≈ 3 at 0–50 m came entirely from each park counting itself. Corrected, there is no other park within 50 m of any sampled park.
- **Commercial land sits close to parks** (F ≈ 3.9 at 0–50 m), and **parks appear around commercial parcels from 50 m outward**. The original denominator understated both.
- **Residential parcels are surrounded by residential land at about the study-area rate** (F ≈ 1). That is expected when about 80% of the non-vacant area is residential.

### What a complete version would need

- Compute enrichment factors per parcel, with the study area as denominator, then map the per-parcel values. This produces a real neighbourhood-effect surface instead of one number per class.
- Use **all parcels**, or a random stratified sample with enough parcels of each class, and report n per curve.
- Assign neighbours by network distance to the parcel itself (for example, to the nearest point on the street that parcel fronts), not via the nearest node, and exclude the central parcel.
- Compute the service areas once per street node and reuse them; the per-parcel Dijkstra runs are what made the full run too slow to finish.

## Repository contents

```
original/land_use_enrichment_factor_2024.ipynb                  # Part 1, as written in 2024 (Persian explanations)
original/part2_network/network_enrichment_factor_2024.ipynb    # Part 2, the exploratory working notebook as written in 2024
corrected/correct_denominator.py                                 # Part 1: exact denominator correction from saved outputs
corrected/part2_network_sample_ef.py                             # Part 2: sample enrichment factors recomputed from saved outputs
figures/                                                         # figures used in this README
```

The land-use rasters and vector data are not included. To reproduce the corrected figures:

```bash
pip install -r requirements.txt
python corrected/correct_denominator.py
python corrected/part2_network_sample_ef.py
```

## References

- Verburg, P. H., de Nijs, T. C. M., Ritsema van Eck, J., Visser, H., & de Jong, K. (2004). A method to analyse neighbourhood characteristics of land use patterns. *Computers, Environment and Urban Systems*, 28(6), 667–690.
- White, R., & Engelen, G. (1993). Cellular automata and fractal urban form: a cellular modelling approach to the evolution of urban land-use patterns. *Environment and Planning A*, 25(8), 1175–1199.

## License

[MIT](LICENSE) © 2024–2026 Kimya Kheirkhah · [github.com/KimyaKheirkhah](https://github.com/KimyaKheirkhah)

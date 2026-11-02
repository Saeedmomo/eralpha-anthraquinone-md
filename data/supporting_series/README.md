# Supporting screening series

These processed records correspond to V9 Section 3.6, main Figure 11 and SI Tables S4–S5 and Figure S4.

| File | Contents |
|---|---|
| `anchor_engagement_V9.csv` | 53 energy records; structural/contact fields are unavailable for 2OUZ and 3ERT. The other 51 complexes support the primary contact analysis. |
| `engagement_definitions_V9.csv` | All eighteen definitions evaluated on the 51 coordinate-supported systems, using Pearson correlation of hydroxyl count. |
| `statistics_V9.json` | Primary, preparation-filtered and energy-stable subsets with 51, 49 and 18 systems. |

Engagement sums the selected aryl-oxygen hydrogen-bond and ionic-contact occupancies. The two-anchor flag requires each marginal residue occupancy to exceed 0.05; it does not establish simultaneous contact in the same frame. The energy-stability filter concerns the energy time series. Prime scores and contact occupancies retain their separate definitions.

The release companion archive contains the original processed source tables and scripts needed to reconstruct these corrected values. Its older records are retained as provenance; the V9 tables and subset definitions are the current interpretation.

# Figure 1 data and panel F update

`Figure1_V11.png` is the revised manuscript figure. Panels A–E, the first-replicate span bars and their within-trajectory standard deviations are unchanged. Only the illustrative reference-width guides in panel F have changed.

The rounded reference width is 14.22 Å, from 14.216 Å in the first 2,6-DHAQ trajectory. At that reference width, the zero reach margin is `s = W - 2a = 8.62 Å`, with `a = 2.8 Å`; the necessary contact-bound reference is `s = W - 2d0 = 7.22 Å`, with `d0 = 3.5 Å`. These guides are examples at one reference width, not universal bridging cutoffs. The full instantaneous geometric condition uses each frame's width and hydroxyl span.

`s4_an_geometry.csv` contains the four original first-replicate geometry rows. `s4_an_hydroxyl.csv` supplies the unchanged chelate data. `msfig1.py` and `msfig_common.py` preserve the original plotting implementation. The local `update_figure1f.py` redirects data paths to this folder and produces PNGs only; it reconstructs the original panel F before applying the two guide annotations to `Figure1_current.png`. `figure1f_validation.json` records unchanged panels, geometry and source hashes.

To regenerate the update:

```bash
python -m pip install matplotlib numpy pillow
python update_figure1f.py
```

Exact pixel reproduction requires the Arial font and the plotting runtime used for the archived figure. On Windows, Arial is supplied by the system; the script also registers Arial from `/mnt/c/Windows/Fonts` when run under WSL. The regression assertion intentionally fails if the original panel cannot be reproduced exactly in a different font or rendering environment. RDKit is not required because the chemical drawings are preserved from the original artwork.

# Hydroxyl availability and geometry figure

`Figure1.png` shows chemical structures, chelation occupancies and first-replicate hydroxyl spans with within-trajectory standard deviations. The reference-width guides in panel F illustrate the reach-margin and instantaneous-contact equations at a cavity width of 14.22 Å. They are not universal cutoffs.

`Figure1_original.png` supplies the unchanged chemical-structure and chelation artwork for panel-specific regeneration. The plotting script reads the accompanying geometry and hydroxyl tables, preserves panels A–E, and verifies the original bar data and geometry. It writes PNG files only.

Run `python update_figure1f.py` with NumPy, Matplotlib, Pillow and Arial available. Original whole-figure generation in `msfig1.py` additionally requires RDKit. The original scientific source and input tables are retained for reproducibility.

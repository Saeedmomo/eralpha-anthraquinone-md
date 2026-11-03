# V11 instantaneous geometric-bound verification

The per-frame table contains 11,022 saved frames from the same twenty-two two-hydroxyl antagonist trajectories used in the manuscript. Frames within trajectories are correlated observations, not independent replicates.

Run the standard-library verification with Python 3:

```text
python verify_v11_frame_bound.py v11_frame_bound_per_frame.csv
```

For each frame, the hard geometric exclusion is `cavity_width_A > instantaneous_hydroxyl_span_A + 7.0`. Direct anchor engagement is `max(GLU353_min_ligand_O_distance_A, HIS524_min_ligand_O_distance_A) < 3.5`. The verifier reproduces 4,552 excluded frames and zero frames with simultaneous engagement. It also checks all twenty-two 501-frame sequences and the agreement of reconstructed instantaneous-span means with the archived run means. No saved anchor distance equals the 3.5 Å cutoff.

The system-mean-span comparison also gives 4,552 excluded frames and zero bridges; this coincidentally equal aggregate count should not be substituted for the instantaneous calculation. The run-mean comparison gives 4,554 and zero. All three are reported in the audit JSON.

For coordinate-level reconstruction, run `frame_bound_audit.py --data-root PATH_TO_MD_S4 --output-dir OUTPUT_DIRECTORY` using NumPy and SciPy. The script reads each run's `an/dry_nofit.nc` and `an/structure.json`, the minimum of the four cavity-width columns in `stageA/span_check/RUN/span.dat`, and the saved `an/triad_GLU353.dat` and `an/triad_HIS524.dat` measurements. One-based oxygen indices in `structure.json` identify the hydroxyl pair. Reconstructed means and SDs match the original records within 0.000001 Å. Inputs and saved trajectories are never modified, and no simulation is rerun.

Coordinate trajectories are required only for reconstructing the derived table from original coordinates. The released per-frame table and standard-library verifier are sufficient to reproduce the reported 4,552/0 count. This comparison concerns sampled geometric engagement and does not estimate affinity or drug efficacy.

SHA-256 of `v11_frame_bound_per_frame.csv`: `f85cbacf14cc542e7dad769f2458e4b1381f29a22632dfdd281af00cbbd6edc9`.

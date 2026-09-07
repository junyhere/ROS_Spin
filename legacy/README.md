# Legacy semiconductor/two-qubit materials

These files are the original historical scripts, retained unchanged alongside
the legacy `dataset/` tree. They model semiconductor spin shuttling, Bell-state
parity, gate noise, CZ depth, transport weighting, or numerical errors in those
calculations. They are not active anthracycline chemistry.

In particular, singlet/triplet two-spin labels cannot be translated into the
doublet/quartet sectors of a spin-1/2 semiquinone plus spin-1 O2 encounter.
Likewise, `f_ow`/`f_fb`, semiconductor damping, Bell parity, and CZ depth cannot
parameterize chemical relaxation, preparation, reaction, or escape. Historical
scripts may require their original dependency versions and are kept for source
inspection, not supported execution.

Use `../spin_chemistry.py`, `../ROS.py`, `../sensitivity_analysis.py`, and
`../paper_analysis.py` for the active six-state sensitivity framework.

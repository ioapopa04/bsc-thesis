# runs/ — the inputs of every simulation, copied from /scratch/ipopa/thesis/<RUN> on sbrinz

One folder per dataset, same layout as on sbrinz, with only what is needed to rebuild and rerun:
`packages/` and `build/parameter/` (ligand parameter packages — reuse them, never re-parameterise),
`prep/` (prepared complexes + atom maps), `*.config` (build-top / scaler / REST2 / parameterize),
`scripts/*.sbatch|py` (build, run, check, analysis), start SDFs, torsion files, job lists.

Not copied (rebuildable or too big): `REST2-run*/` (trajectories), `built.*`, `min/`, `input/`,
analysis outputs, logs, `test/` runs, `_superseded*/`.

To use one on a new machine: `cp -a runs/<RUN> $THESIS/<RUN>`, rewrite paths (../SETUP.md §2),
then rebuild with its `scripts/build.sbatch` and run with `scripts/rest2.sbatch`.
Which run is which: ../SETUP.md §5.  folic_acid_mace_omol/ = Shu-Yu's bespoke MACE-OMOL System XML.

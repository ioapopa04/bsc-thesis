Bachelor thesis -- everything in one place (moved here 2026-10-15)

bsc-thesis/                git repo (github.com:ioapopa04/bsc-thesis): scripts, results, crystal ligands
FOLn-REST2/                folate in water, first REST2 run
FOLn-REST2-50ns/           folate in water, REST2 100 ns, 8 replicas (original REST2)
FOL-complex-REST2/         folate bound to 4LRH (FRalpha) and 3QL3 (DHFR), original REST2
FOL-REST2-v4/              folate water / 4LRH / 3QL3 with REST2 v4 (COO-, amide heated)
REST2-PARA/                paracetamol in water, original REST2, 10 ns
PARA-REST2-v4/             paracetamol in water, REST2 v4, 10 ns
REST2-conformers/          5 small molecules (4ibi 4ibk 4j7d 5mwh 5ov8) in water
REST2-conformers-complex/  the same 5 bound to their proteins, 40 ns
tmp-rest2-test/            early REST2 test
folenv.sh                  shell setup for the FOL workflow

Left OUTSIDE on purpose (absolute paths baked in): /scratch/ipopa/software (conda envs, MD-tools
checkouts) and /scratch/ipopa/md-data (MD-tools data, MD_DATA).
Run records (restart.json, *.log, *.out) still show the old paths; they are history, not inputs.

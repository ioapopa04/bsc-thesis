---
name: thesis-folder
description: All bachelor-thesis folders were moved under /scratch/ipopa/thesis/ on 2026-10-15; software/ and md-data/ stay at /scratch/ipopa
metadata:
  type: project
---
User asked to gather everything thesis-related in one parent folder (other projects clutter /scratch/ipopa).
Moved: bsc-thesis repo, FOLn-REST2(-50ns), FOL-complex-REST2, FOL-REST2-v4, REST2-PARA, PARA-REST2-v4, REST2-conformers(-complex), tmp-rest2-test, folenv.sh -> /scratch/ipopa/thesis/.
Not moved: /scratch/ipopa/software (conda envs have baked-in prefixes), /scratch/ipopa/md-data.
**How to apply:** new thesis runs go in /scratch/ipopa/thesis/. Script paths were rewritten; run records (restart.json, logs) and REST2-conformers/extend_all.sh (user: don't touch) still have old paths.

#!/bin/bash
# Submit the 12 scenarios: 3 life histories x {no insurance, insurance} x {no relief, relief}
# on the CQLS HPC. Run from anywhere on the cluster:  bash run_scenarios.sh
# Results: results/V2/scenarios/<LH>_<ins|noins>_<rel|norel>/
# Logs:    scenario-<LH>-ins<0|1>-rel<0|1>-StdOut/  (in the directory you run this from)

DIRPATH=/home/ceoas/bucknejo/github/
NODE=amaterasu01
PROCS=20

# grids (same for every scenario so the results are comparable)
NKAPPA=20
NETA=20
NS=20
NB=20
NDELTAB=15
NQUAD=10
MAXITER=100

for lh in RA_short RA_intermidiate RA_long; do
  for ins in 0 1; do
    for rel in 0 1; do
      hqsub -P $PROCS "julia --threads=$PROCS ${DIRPATH}FisheriesInsuranceModel/src/run_scenario.jl $lh $ins $rel $NKAPPA $NETA $NS $NB $NDELTAB $NQUAD $MAXITER" \
        -r scenario-$lh-ins$ins-rel$rel-StdOut -q ceoas@$NODE
    done
  done
done

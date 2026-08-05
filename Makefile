################################################
### Runs batched analyysis of the fishiersy  ###
### insurnace decison model on the CQLS HPC  ###
################################################
#
#
# Must end in a slash: every target below concatenates it directly, as in
# $(DIRPATH)FisheriesInsuranceModel/src/run.jl
DIRPATH=/home/ceoas/bucknejo/github/
#
#
PROCS=20
#
#
NODE=amaterasu01
#
#
PARAMS=RA_intermidiate
#
#
NKAPPA=15
#
#
NETA=10
#
#
NS=20
#
#
NB=20
#
#
NDELTAB=15
#
#
NQUAD=10
#
#
all:
	@echo ""
	@echo "    This folder runs the fishiersy insurance model on the CQLS HPC "
	@echo ""
	@echo ""
	@echo "    Please type.....   "
	@echo ""
	@echo "    	make run "
	@echo "	make run_params  PARAMS=<PARAMS> "
	@echo "	make clean  PARAMS=<PARAMS> "
	@echo ""
	@echo "    index insurance model (run on the cluster, in this order): "
	@echo ""
	@echo "	make setup_index                        instantiate the env, once "
	@echo "	make index_jobs                         list what would be submitted "
	@echo "	make run_index_params PARAMS=<SCENARIO>  one job, for a smoke test "
	@echo "	make run_index_no_insurance             2 jobs "
	@echo "	make run_index                          14 jobs "
	@echo "	make sync_index                         cluster -> local (run locally) "
#
#
run:
	@for params in `ls -1 $(DIRPATH)FisheriesInsuranceModel/parameters | grep .jl |  sed 's/.jl//g'`; do \
    hqsub -P $(PROCS) "julia --threads=$(PROCS) $(DIRPATH)FisheriesInsuranceModel/src/run.jl  $${params} $(NKAPPA) $(NETA) $(NS) $(NB) $(NDELTAB) $(NQUAD)" -r solve-$${params}-StdOut -q ceoas@$(NODE);\
	hqsub -P $(PROCS) "julia --threads=$(PROCS) $(DIRPATH)FisheriesInsuranceModel/src/run_no_insurance.jl  $${params}" -r solve_no_insurance-$${params}-StdOut -q ceoas@$(NODE);\
    done;
#
#
run_params:
	hqsub -P $(PROCS) "julia --threads=$(PROCS) $(DIRPATH)FisheriesInsuranceModel/src/run.jl  ${PARAMS} $(NKAPPA) $(NETA) $(NS) $(NB) $(NDELTAB) $(NQUAD)" -r solve-${PARAMS}-StdOut -q ceoas@$(NODE);
	hqsub -P $(PROCS) "julia --threads=$(PROCS) $(DIRPATH)FisheriesInsuranceModel/src/run_no_insurance.jl  ${PARAMS} $(NKAPPA) $(NETA) $(NS) $(NB) $(NDELTAB) $(NQUAD)" -r solve_no_insurance-${PARAMS}-StdOut -q ceoas@$(NODE)
#
#
clean:
	rm -r -f *-StdOut


################################################
### Index insurance model  (index_model/)
###
### Kept separate from the base model targets above, which are left
### untouched.  Three differences from the base model contract:
###
###   1. the runners take NO grid size arguments - sizes come from
###      index_model/config.toml, so every job in a sweep is on the
###      same grid by construction (NKAPPA ... NQUAD above apply to
###      the base model targets only);
###   2. the no insurance solve is submitted once per BASE parameter
###      set, not once per scenario (see IDX_NO_INSURANCE below);
###   3. paths are named once, in IDX_ROOT, rather than repeated.
################################################

# The checkout on the cluster: /home/ceoas/bucknejo/github/FisheriesInsuranceModel
IDX_ROOT=$(DIRPATH)FisheriesInsuranceModel
IDX_PARAMS=$(IDX_ROOT)/index_model/parameters
IDX_HOST=bucknejo@hpc.ceoas.oregonstate.edu

# Scenario names, read off the parameter directory.  generate_scenarios.jl
# writes those files and is not a scenario itself, so it must not be submitted.
IDX_SCENARIOS=$$(ls -1 $(IDX_PARAMS) | grep '\.jl$$' | grep -v -e generate_scenarios -e universal | sed 's/\.jl$$//')

# One no insurance solve per base parameter set.  The no insurance value
# function is invariant to (pr, re) and (cf, cv), so all seven scenarios on a
# base share ONE output directory - submitting one job per scenario would race
# seven writers onto the same sol.jld2.  Any scenario on the base works as the
# argument; the runner resolves it to no_insurance_<base>/.
IDX_NO_INSURANCE=no_RA_int_pr100_re100 RA_int_pr100_re100


# Run once, on the head node, before submitting anything.  Kept out of the
# workers deliberately: dozens of jobs racing on the depot to fetch the
# unregistered ValueFunctionIterations.jl git dependency is a real failure mode.
setup_index:
	julia --project=$(IDX_ROOT) -e 'using Pkg; Pkg.instantiate(); Pkg.status()'

# Print what run_index would submit, without submitting it.  Cheap check that
# the scenario list and the cluster paths are what you think they are.
index_jobs:
	@echo "root       $(IDX_ROOT)"
	@ls $(IDX_ROOT)/index_model/config.toml >/dev/null && echo "config     ok"
	@echo "no insur.  $(IDX_NO_INSURANCE)"
	@echo "scenarios:"
	@for params in $(IDX_SCENARIOS); do echo "    $${params}"; done
	@echo "total      $$(echo $(IDX_SCENARIOS) | wc -w) insurance + $$(echo $(IDX_NO_INSURANCE) | wc -w) no insurance jobs"

# 2 jobs.  Independent of run_index, so the order does not matter.
run_index_no_insurance:
	@for base in $(IDX_NO_INSURANCE); do \
	  hqsub -P $(PROCS) "julia --project=$(IDX_ROOT) --threads=$(PROCS) $(IDX_ROOT)/index_model/run_no_insurance.jl $${base}" -r solve-index-no-insurance-$${base}-StdOut -q ceoas@$(NODE); \
	done;

# 14 jobs, one per scenario cell.
run_index:
	@for params in $(IDX_SCENARIOS); do \
	  hqsub -P $(PROCS) "julia --project=$(IDX_ROOT) --threads=$(PROCS) $(IDX_ROOT)/index_model/run.jl $${params}" -r solve-index-$${params}-StdOut -q ceoas@$(NODE); \
	done;

# One scenario, for the smoke test: submit this alone and check the wall clock
# and the output directory before committing to the full sweep.
run_index_params:
	hqsub -P $(PROCS) "julia --project=$(IDX_ROOT) --threads=$(PROCS) $(IDX_ROOT)/index_model/run.jl ${PARAMS}" -r solve-index-${PARAMS}-StdOut -q ceoas@$(NODE)

run_index_no_insurance_params:
	hqsub -P $(PROCS) "julia --project=$(IDX_ROOT) --threads=$(PROCS) $(IDX_ROOT)/index_model/run_no_insurance.jl ${PARAMS}" -r solve-index-no-insurance-${PARAMS}-StdOut -q ceoas@$(NODE)

# The whole sweep, no insurance first so a mistake in the cheap jobs surfaces
# before the expensive ones queue.
run_index_model: run_index_no_insurance run_index

# Cluster -> local, one directional.  results_index/ is gitignored.
sync_index:
	rsync -avz --partial --info=progress2 $(IDX_HOST):$(IDX_ROOT)/results_index/ ./results_index/
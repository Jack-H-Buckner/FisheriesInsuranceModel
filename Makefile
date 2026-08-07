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
	@echo ""
	@echo "	make index_status                       one line per job, from its log "
	@echo "	make index_log PARAMS=<SCENARIO>        full log for one job "
	@echo "	make index_errors                       only the jobs that failed "
	@echo "	make index_results                      what has landed in results_index/ "
	@echo ""
	@echo "    index insurance analysis - bankruptcy rates, on the cluster: "
	@echo ""
	@echo "	make bankruptcy_jobs                    list what would be submitted "
	@echo "	make run_index_bankruptcy_quick         all scenarios, tiny, ~2 min "
	@echo "	make run_index_bankruptcy               one job per scenario "
	@echo "	make bankruptcy_status                  one line per job, from its log "
	@echo "	make bankruptcy_errors                  only the jobs that failed "
	@echo ""
	@echo "    ... then, locally: "
	@echo ""
	@echo "	make sync_index_analysis                cluster -> local "
	@echo "	make merge_index_bankruptcy             parts -> index_bankruptcy.csv "
	@echo ""
	@echo "	make analysis_index_bankruptcy          the same run, but all in one "
	@echo "	                                        process on this machine "
	@echo "	make test_index                         the index_model test suite "
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
	rsync -avz $(IDX_HOST):$(IDX_ROOT)/results_index/ ./results_index/


################################################
### Index insurance analysis: bankruptcy rates (demand metric 3)
###
### Submitted on the cluster, like the solves, but with a different
### shape.  The solves are one long job per scenario; this is one SHORT
### job per scenario, and the parallel unit inside a job is the N
### simulated trajectories rather than the state grid.
###
### It reads results_index/ and never writes it, so it runs on the
### cluster directly after run_index finishes - there is no need to
### sync the 300 MB of solutions down first.  Each job writes its own
### part table; merge_index_bankruptcy combines them.
###
###   on the cluster        make run_index_bankruptcy
###                         make bankruptcy_status
###   locally               make sync_index_analysis
###                         make merge_index_bankruptcy
###
### A scenario whose solve has not landed yet is reported and skipped
### by its job, so a partial sweep still produces a partial table.
################################################

IDX_MANUSCRIPT=$(IDX_ROOT)/manuscript_index
IDX_PARTS=$(IDX_MANUSCRIPT)/parts

# Scenarios to analyse.  The plan computes the demand metrics on the risk
# neutral sets, so this is the scenario list filtered to no_RA_*; drop the
# `grep '^no_RA'` to include the RA_* sets as well.
IDX_ANALYSIS_SCENARIOS=$$(ls -1 $(IDX_PARAMS) | grep '\.jl$$' | grep -v -e generate_scenarios -e universal | grep '^no_RA' | sed 's/\.jl$$//')

# Extra arguments for run_bankruptcy.jl, e.g.
#   make run_index_bankruptcy BANKRUPTCY_ARGS="--N=120 --B=4000"
# Run `julia index_model/analysis/run_bankruptcy.jl --help` style options are
# documented in the header of that file.
BANKRUPTCY_ARGS=

# Print what run_index_bankruptcy would submit, without submitting it.
bankruptcy_jobs:
	@echo "root       $(IDX_ROOT)"
	@echo "parts      $(IDX_PARTS)"
	@echo "args       $(BANKRUPTCY_ARGS)"
	@echo "scenarios:"
	@for params in $(IDX_ANALYSIS_SCENARIOS); do \
	  if [ -f $(IDX_ROOT)/results_index/$${params}/sol.jld2 ]; then s=solved; else s="NOT SOLVED"; fi; \
	  printf "    %-34s %s\n" "$${params}" "$${s}"; \
	done
	@echo "total      $$(echo $(IDX_ANALYSIS_SCENARIOS) | wc -w) jobs"

# One job per scenario.  Each is threaded over its N trajectories, so PROCS
# threads is the right request even though the job is short.
run_index_bankruptcy:
	@mkdir -p $(IDX_PARTS)
	@for params in $(IDX_ANALYSIS_SCENARIOS); do \
	  hqsub -P $(PROCS) "julia --project=$(IDX_ROOT) --threads=$(PROCS) $(IDX_ROOT)/index_model/analysis/run_bankruptcy.jl $${params} --out=$(IDX_PARTS)/index_bankruptcy_$${params}.csv $(BANKRUPTCY_ARGS)" -r bankruptcy-$${params}-StdOut -q ceoas@$(NODE); \
	done;

# The whole sweep at a resolution that finishes in a couple of minutes.  Submit
# this first: it exercises every scenario's solution end to end, so a missing or
# mismatched solve surfaces before the real run queues.  The numbers it writes
# are a plumbing check, NOT results - overwrite them with the real run.
run_index_bankruptcy_quick:
	@$(MAKE) run_index_bankruptcy BANKRUPTCY_ARGS="--quick $(BANKRUPTCY_ARGS)"

# One scenario:  make run_index_bankruptcy_params PARAMS=no_RA_int_pr60_re100
run_index_bankruptcy_params:
	@mkdir -p $(IDX_PARTS)
	hqsub -P $(PROCS) "julia --project=$(IDX_ROOT) --threads=$(PROCS) $(IDX_ROOT)/index_model/analysis/run_bankruptcy.jl ${PARAMS} --out=$(IDX_PARTS)/index_bankruptcy_${PARAMS}.csv $(BANKRUPTCY_ARGS)" -r bankruptcy-${PARAMS}-StdOut -q ceoas@$(NODE)

# One line per job.  The markers come from run_bankruptcy.jl's own output:
#   "simulating"  -> solutions loaded, trajectories running
#   "wrote"       -> the part table is on disk, job complete
bankruptcy_status:
	@found=0; \
	for d in bankruptcy-*-StdOut; do \
	  [ -d "$$d" ] || continue; found=1; \
	  name=$${d#bankruptcy-}; name=$${name%-StdOut}; \
	  if grep -qs "^ERROR" $$d/* 2>/dev/null; then state=FAILED; \
	  elif grep -qs "^wrote " $$d/* 2>/dev/null; then state=done; \
	  elif grep -qs "skipping" $$d/* 2>/dev/null; then state=SKIPPED; \
	  elif grep -qs "simulating" $$d/* 2>/dev/null; then state=running; \
	  else state=starting; fi; \
	  last=`cat $$d/* 2>/dev/null | tr '\r' '\n' | grep -v '^[[:space:]]*$$' | tail -1 | cut -c1-90`; \
	  printf "%-9s %-34s %s\n" "$$state" "$$name" "$$last"; \
	done; \
	[ $$found = 1 ] || echo "no bankruptcy-*-StdOut directories here - run this where you ran run_index_bankruptcy"

# Full log for one job:  make bankruptcy_log PARAMS=no_RA_int_pr60_re100
bankruptcy_log:
	@cat bankruptcy-$(PARAMS)-StdOut/* 2>/dev/null | tr '\r' '\n' | tail -n 60 || \
	  echo "no log for $(PARAMS)"

bankruptcy_errors:
	@for d in bankruptcy-*-StdOut; do \
	  [ -d "$$d" ] || continue; \
	  grep -qs "^ERROR" $$d/* 2>/dev/null || continue; \
	  echo "===== $$d"; \
	  grep -A5 "^ERROR" $$d/* 2>/dev/null | head -20; \
	done

# What has landed in manuscript_index/parts/ on the cluster.
bankruptcy_results:
	@ls -1 $(IDX_PARTS) 2>/dev/null | sed 's/^/    /' || echo "nothing in $(IDX_PARTS) yet"


################################################
### Index insurance analysis: local side
################################################

# Cluster -> local, one directional, and small: these are CSVs, not solutions.
sync_index_analysis:
	rsync -avz $(IDX_HOST):$(IDX_MANUSCRIPT)/ ./manuscript_index/

# parts/*.csv -> manuscript_index/index_bankruptcy.csv, in sweep order.
# Cheap and idempotent, so it is fine to rerun while jobs are still landing.
merge_index_bankruptcy:
	julia --project=. index_model/analysis/merge_bankruptcy.jl

# IDX_THREADS is the parallel width of the bankruptcy simulation when it is run
# here rather than submitted; set it to the number of cores on this machine.
IDX_THREADS=8

# The whole sweep in ONE process, no scheduler, writing the merged CSV directly.
# Use this locally after `make sync_index`; on the cluster prefer the fan out
# above, which is the same calculation spread over one job per scenario.
analysis_index_bankruptcy:
	julia --project=. --threads=$(IDX_THREADS) index_model/analysis/run_bankruptcy.jl $(BANKRUPTCY_ARGS)

analysis_index_bankruptcy_quick:
	julia --project=. --threads=$(IDX_THREADS) index_model/analysis/run_bankruptcy.jl --quick $(BANKRUPTCY_ARGS)

analysis_index_bankruptcy_params:
	julia --project=. --threads=$(IDX_THREADS) index_model/analysis/run_bankruptcy.jl ${PARAMS} $(BANKRUPTCY_ARGS)

test_index:
	julia --project=. --threads=$(IDX_THREADS) index_model/test/test_transitions.jl
	julia --project=. --threads=$(IDX_THREADS) index_model/test/test_grids.jl
	julia --project=. --threads=$(IDX_THREADS) index_model/test/test_runners.jl
	julia --project=. --threads=$(IDX_THREADS) index_model/test/test_analysis.jl
	julia --project=. --threads=$(IDX_THREADS) index_model/test/test_coverage.jl
	julia --project=. --threads=$(IDX_THREADS) index_model/test/test_bankruptcy.jl


################################################
### Job logs
###
### hqsub collects each job's output in a <name>-StdOut DIRECTORY,
### created in the directory make was run from - so run these from the
### same place you ran run_index (and `make clean` deletes them).
### The runners print a fixed set of markers, which is what
### index_status greps for:
###
###   "solving (maxiter = " -> the grid is built, VFI has started
###   "solved in"           -> VFI finished, writing the solution
###   "wrote diagnostic"    -> everything on disk, job complete
################################################

# One line per job: state, scenario, and the last thing it printed (the
# ProgressMeter bar is \r separated, so it is split back into lines to show the
# current percentage and ETA rather than one unreadable blob).
index_status:
	@found=0; \
	for d in solve-index-*-StdOut; do \
	  [ -d "$$d" ] || continue; found=1; \
	  name=$${d#solve-index-}; name=$${name%-StdOut}; \
	  if grep -qs "^ERROR" $$d/* 2>/dev/null; then state=FAILED; \
	  elif grep -qs "wrote diagnostic plots" $$d/* 2>/dev/null; then state=done; \
	  elif grep -qs "solved in" $$d/* 2>/dev/null; then state=saving; \
	  elif grep -qs "solving (maxiter" $$d/* 2>/dev/null; then state=solving; \
	  else state=starting; fi; \
	  last=`cat $$d/* 2>/dev/null | tr '\r' '\n' | grep -v '^[[:space:]]*$$' | tail -1 | cut -c1-90`; \
	  printf "%-9s %-30s %s\n" "$$state" "$$name" "$$last"; \
	done; \
	[ $$found = 1 ] || echo "no solve-index-*-StdOut directories here - run this where you ran run_index"

# Full log for one job:  make index_log PARAMS=no_RA_int_pr60_re100
index_log:
	@cat solve-index-$(PARAMS)-StdOut/* 2>/dev/null | tr '\r' '\n' | tail -n 60 || \
	  echo "no log for $(PARAMS)"

# Just the failures, with the error text.
index_errors:
	@for d in solve-index-*-StdOut; do \
	  [ -d "$$d" ] || continue; \
	  grep -qs "^ERROR" $$d/* 2>/dev/null || continue; \
	  echo "===== $$d"; \
	  grep -A5 "^ERROR" $$d/* 2>/dev/null | head -20; \
	done

# What has actually landed in results_index/ on the cluster.
index_results:
	@for d in `ls -1 $(IDX_ROOT)/results_index 2>/dev/null`; do \
	  printf "%-34s %s\n" "$$d" "`ls $(IDX_ROOT)/results_index/$$d | tr '\n' ' '`"; \
	done
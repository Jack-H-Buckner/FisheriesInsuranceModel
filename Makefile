################################################
### Runs batched analyysis of the fishiersy  ###
### insurnace decison model on the CQLS HPC  ###
################################################
#
#
DIRPATH=/home/ceoas/bucknejo/
#
#
PROCS=40
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


#
#
# The index runners take no size arguments: grid sizes come from
# index_model/config.toml, so every job in the sweep uses the same grid.
# NKAPPA ... NQUAD above apply to the base model targets only.
run_index_model:
	@for params in `ls -1 $(DIRPATH)FisheriesInsuranceModel/index_model/parameters | grep .jl |  sed 's/.jl//g'`; do \
    hqsub -P $(PROCS) "julia --threads=$(PROCS) $(DIRPATH)FisheriesInsuranceModel/index_model/run.jl  $${params}" -r solve-$${params}-StdOut -q ceoas@$(NODE);\
	hqsub -P $(PROCS) "julia --threads=$(PROCS) $(DIRPATH)FisheriesInsuranceModel/index_model/run_no_insurance.jl  $${params}" -r solve_no_insurance-$${params}-StdOut -q ceoas@$(NODE);\
    done;
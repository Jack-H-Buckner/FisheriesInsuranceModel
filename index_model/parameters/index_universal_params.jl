#####################################################################
#####################################################################
### Shared contract settings for the index insurance scenarios.
###
### Event frequency (p_low) and duration (p_stay) are deliberately NOT
### set here.  set_index_T reads them out of the base parameter set's
### 2x2 T, so layering the same scenario grid on a different base file
### (e.g. no_RA_intermidiate_higher_frequency.jl) needs no changes here.
### Jack H. Buckner, Oregon State University, 2026
#####################################################################
#####################################################################

index_cf = 1e-6 # fixed insurance cost; base model default (the notebook used 0.0)
index_cv = 0.0  # variable cost / load; base model default (the notebook used 0.1)

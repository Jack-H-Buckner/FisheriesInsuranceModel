include("../src/derived_parameters.jl")

# constants to facilitate computation


function joint_transition_from_matrix(T_mx0,T_mx1,T_x)
    E0 = [1.0 0.0; 0.0 0.0]
    E1 = [0.0 0.0; 0.0 1.0]
    return kron(T_mx0, E0*T_x) .+ kron(T_mx1, E1*T_x)
end


function joint_transition(p_low,pr, re, p_stay)
    p_index = re*p_low/pr 
    b0      = p_low*(1-re)/(1-p_index) 

    T_mx0 = [(1-b0) (1-p_stay); b0 p_stay]
    T_mx1 = [(1-pr) (1-p_stay); pr p_stay] 
    T_x   = [(1-p_index) (1-p_index);  p_index p_index]

    joint_transition_from_matrix(T_mx0,T_mx1,T_x)
end 
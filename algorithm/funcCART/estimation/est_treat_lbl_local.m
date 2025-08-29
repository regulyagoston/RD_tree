%% Estimates the treatment effect using local polynomial estimator
%
% Input:
%   node - nodeProp output
%       contains the node/leaf which needs to be estimated
%   data - sample outect
%       contains all the data for the estimation
%   optTree - optCART output
%       contains relevant information about what to estimate
%   tr - 1x1 logical value
%       true  - it uses the training sample
%       false - it uses the estiation sample
%
% Output
%   out - structure
%       output from rdrobust
%       p_p_j      - share of treated units in the leaf
%       p_m_j      - share of untreated units in the leaf

function out = est_treat_lbl_local( node , data , optTree , tr )


%% Get the data to work with (training or estimation sample)
if tr
    [ Y , X , W , ~ , cl ] = get_node_tr_Sample( node , data );
else
    [ Y , X , W , ~ , cl ] = get_node_est_Sample( node , data );
end

if strcmp( optTree.type , 'CLATE' )
    fuzzy = W;
else
    fuzzy = [];
end
if isnan( cl )
    cl = [];
end

%% Local polinomial estimation
if ( strcmp( optTree.bw_type, 'leaf-by-leaf' ) )
    bw_j = rdbwselect( Y , X( : , 1 ) , 'c', data.c, 'fuzzy', fuzzy, 'p', optTree.orderPolinomial, ...
        'vce', optTree.varEstimator, 'bwselect', optTree.bwselect );
else
    bw_j = struct;
    bw_j.h_l = optTree.bw;
    bw_j.h_r = optTree.bw;
    bw_j.b_l = optTree.bw_b;
    bw_j.b_r = optTree.bw_b;
end

out_raw = rdrobust( Y , X( : , 1 ) , data.c , bw_j, fuzzy , cl , optTree.orderPolinomial,...
                optTree.varEstimator, optTree );

%% Reconstruct output for ML
n_ell = length( Y );
out.tau_i = repmat( out_raw.tau_bc, [ n_ell, 1 ] );
out.tau_j = out_raw.tau_bc;
out.tau_raw_j = out_raw.tau_cl;
out.var_tau_cl_j = out_raw.tau_V_cl;
out.var_tau_bc_j = out_raw.tau_V_rb;
out.tau_se_j = out_raw.tau_se_rb;
out.bias_j = out_raw.bias_r - out_raw.bias_l;
% Note here bandiwdth assumed to be the same below and above
out.h_j = bw_j.h_l;
out.v_j = out_raw.tau_V_cl*n_ell*out.h_j;
out.bw_j = bw_j;

    
%% Add further information            
% Honest and treatment type values
if optTree.honest && any( strcmp( optTree.type , {'CATE','CLATE'} ) )
    out.n_p_j_tr  = sum( node.logID_tr &  data.ITT_tr );
    out.n_m_j_tr  = sum( node.logID_tr & ~data.ITT_tr );
    out.n_p_j_est = sum( node.logID_est &  data.ITT_est );
    out.n_m_j_est = sum( node.logID_est & ~data.ITT_est );
    out.p_p_j_tr  = out.n_p_j_tr ./ ( out.n_p_j_tr + out.n_m_j_tr );
    out.p_m_j_tr  = 1 - out.p_p_j_tr;
    out.p_p_j_est = out.n_p_j_est ./ ( out.n_p_j_est + out.n_m_j_est );
    % If no observation in estimation sample, assume it has the same
    % percentage as in the training sample
    if isnan( out.p_p_j_est )
        out.p_p_j_est = out.p_p_j_tr;
    end
    out.p_m_j_est = 1 - out.p_p_j_est;
    if strcmp( optTree.type, 'CLATE' )
        out.tau_y = out_raw.tau_Y_bc;
        %out.tau_y_se = tau_V_Y_rb;
        out.tau_t = out_raw.tau_T_bc;
        out.tau_t_se = sqrt(out_raw.tau_V_T_rb);
    end
end           

end


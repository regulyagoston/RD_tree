%% Treatment criterion using honest EMSE criterion from Athey-Imbens (2016)

function sQ = treat_mean_lbl_ai_crit( leaf_nodes , obj_sample , obj_optCART )

[ tau_i , sigma2_p_j , sigma2_m_j ] = collectInputs( leaf_nodes , true , 'tau_i' , 'sigma2_p_j' , 'sigma2_m_j' );

% Additional information to be used
N_tr  = obj_sample.n_tr;
N_est = obj_sample.n_est;
p     = obj_sample.p_tr;

% Criterion
sQ = - obj_optCART.criteria_weight .* ...
         1./ N_tr .* sum( tau_i .^ 2 ) ...
     + ( 1 - obj_optCART.criteria_weight ) .* ...
       ( 1 ./ N_tr + 1 ./ N_est)   .* sum( sigma2_p_j ./ p + sigma2_m_j ./ ( 1 - p ) ); 

end
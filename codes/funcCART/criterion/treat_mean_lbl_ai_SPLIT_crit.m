%% Treatment criterion using honest EMSE criterion from Athey-Imbens (2016)

function sQ = treat_mean_lbl_ai_SPLIT_crit( split_nodes , obj_sample , obj_optCART )

% Additional information to be used
p     = obj_sample.p_tr;

% Get the estimated values for the treatment and variances
tau_i      = [ split_nodes( 1 ).est_struct.tau_i ; split_nodes( 2 ).est_struct.tau_i ];
n_tr       = numel( tau_i );
sigma2_p_j = [ split_nodes( 1 ).est_struct.sigma2_p_j , split_nodes( 2 ).est_struct.sigma2_p_j ];
sigma2_m_j = [ split_nodes( 1 ).est_struct.sigma2_m_j , split_nodes( 2 ).est_struct.sigma2_m_j ];

% Criterion
sQ = - obj_optCART.criteria_weight .* ...
         1 ./ n_tr .* sum( tau_i .^ 2 ) ...
     + ( 1 - obj_optCART.criteria_weight ) .* ( 1 + obj_sample.train_to_est_ratio ) .* ...
         1 ./ n_tr .* sum( sigma2_p_j ./ p + sigma2_m_j ./ ( 1 - p ) ); 

end
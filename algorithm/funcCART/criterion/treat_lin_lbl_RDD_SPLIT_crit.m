%% Treatment criterion using honest EMSE criterion tailored for RDD
%
% N_tr  - scalar
%   number of all obs in training sample
% N_est - scalar
%   number of all obs in estimation sample
% tau_i - N_tr x 1 vector
%   treatment for individual i
% sigma2_p_j - 1 x 1 x #\Pi vector
%   V[Y|X,Z in ell] of treated for leaf j
% sigma2_m_j - 1 x 1 x #\Pi vector
%   V[Y|X,Z in ell] of untreated for leaf j
% M_p_j - p+1 x p+1 x #\Pi matrix
%   Var-Covar matrix for covariates of treated for leaf j
% M_m_j - p+1 x p+1 x #\Pi matrix
%   Var-Covar matrix for covariates of untreated for leaf j
% p_p_j - 1 x 1 x #\Pi vector
%   leaf share of the treated units
% p_m_j - 1 x 1 x #\Pi vector
%   leaf share of the untreated units

function sQ = treat_lin_lbl_RDD_SPLIT_crit( split_nodes , obj_sample , obj_optCART )

%error();

% Additional information to be used
N_tr  = obj_sample.n_tr;
N_est = obj_sample.n_est;

% Get the estimated values for the treatment and variances
tau_i      = [ split_nodes( 1 ).est_struct.tau_i ; split_nodes( 2 ).est_struct.tau_i ];
sigma2_p_j = cat( 3 , split_nodes( 1 ).est_struct.sigma2_p_j , split_nodes( 2 ).est_struct.sigma2_p_j );
sigma2_m_j = cat( 3 , split_nodes( 1 ).est_struct.sigma2_m_j , split_nodes( 2 ).est_struct.sigma2_m_j );
Mi_p_j     = cat( 3 , split_nodes( 1 ).est_struct.Mi_p_j , split_nodes( 2 ).est_struct.Mi_p_j );
Mi_m_j     = cat( 3 , split_nodes( 1 ).est_struct.Mi_m_j , split_nodes( 2 ).est_struct.Mi_m_j );
p_p_j_est  = cat( 3 , split_nodes( 1 ).est_struct.p_p_j_est , split_nodes( 2 ).est_struct.p_p_j_est );
p_m_j_est  = cat( 3 , split_nodes( 1 ).est_struct.p_m_j_est , split_nodes( 2 ).est_struct.p_m_j_est );

         
if strcmp( obj_optCART.varEstimator , 'simple' )
    % Homoscedastic case: var = Sigma * Mi / p
    var_combined = sigma2_p_j .* Mi_p_j ./ p_p_j_est ...
                 + sigma2_m_j .* Mi_m_j ./ p_m_j_est;
elseif any( strcmp( obj_optCART.varEstimator , {'hce-0','hce-1','clustered' } ) )
    % In general cases: var = Mi * Sigma * Mi / p
    var_combined = zeros( size( Mi_p_j , 1 ) , size( Mi_p_j , 2 ) );
    for i = 1 : 2 
        var_combined = var_combined + ...
                       ( Mi_p_j( : , : , i ) * sigma2_p_j( : , : , i ) * Mi_p_j( : , : , i ) ) ./ p_p_j_est( : , : , i ) ...
                     + ( Mi_m_j( : , : , i ) * sigma2_m_j( : , : , i ) * Mi_m_j( : , : , i ) ) ./ p_m_j_est( : , : , i );
    end
else
    error('No such variance implemented for the criterion!')
end       

% Criterion
sQ = -obj_optCART.criteria_weight .* ...
         1./ N_tr .* sum( tau_i .^ 2 ) ...
     + ( 1 - obj_optCART.criteria_weight ) .* ...
       ( ( 1 ./ N_tr + 1./ N_est ) * sum( var_combined( 1 , 1 , : ) ) );

end
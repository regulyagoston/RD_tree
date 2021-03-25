%% Treatment criterion using honest EMSE criterion tailored for sharp RDD
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

function sQ = treat_lin_lbl_RDD_crit( leaf_nodes , obj_sample , obj_optCART )

[ tau_i , sigma2_p_j , sigma2_m_j , Mi_p_j , Mi_m_j ,...
  p_p_j_est , p_m_j_est ] = collectInputs( leaf_nodes , true , ...
                                                        'tau_i' , 'sigma2_p_j' ,...
                                                        'sigma2_m_j' , 'Mi_p_j' , 'Mi_m_j' ,...
                                                        'p_p_j_est' , 'p_m_j_est' );

% Additional information to be used
N_tr  = obj_sample.n_tr;
N_est = obj_sample.n_est;
nPi   = numel( leaf_nodes );
%
l1 = strcmp( obj_optCART.varEstimator , 'simple' ) && ~isprop( obj_sample , 'cl_tr' );

if l1
    % Homoscedastic case: var = Sigma * Mi / p
    var_combined = reshape( sigma2_p_j , [ 1 , 1 , nPi ] ) .* Mi_p_j ./ reshape( p_p_j_est , [ 1 , 1 , nPi ] ) ...
                 + reshape( sigma2_m_j , [ 1 , 1 , nPi ] ) .* Mi_m_j ./ reshape( p_m_j_est , [ 1 , 1 , nPi ] );
elseif any( strcmp( obj_optCART.varEstimator , {'hce-0','hce-1' } ) ) || ~l1
    % In general cases: var = Mi * Sigma * Mi / p
    if nPi > 1
        var_combined = zeros( size( sigma2_p_j( : , : , 1 ) ) );
        for j = 1 : nPi
            var_combined = var_combined ...
                         + ( Mi_p_j( : , : , j ) * sigma2_p_j( : , : , j ) * Mi_p_j( : , : , j ) ) ./ p_p_j_est( j ) ...
                         + ( Mi_m_j( : , : , j ) * sigma2_m_j( : , : , j ) * Mi_m_j( : , : , j ) ) ./ p_m_j_est( j );    
        end
    else
        var_combined = ( Mi_p_j * sigma2_p_j * Mi_p_j ) ./ p_p_j_est ...
                     + ( Mi_m_j * sigma2_m_j * Mi_m_j ) ./ p_m_j_est;    
    end
else
    error('No such variance implemented for the criterion!')
end
% Criterion
sQ = - obj_optCART.criteria_weight .* ...
         1./ N_tr .* sum( tau_i .^ 2 ) ...
     + ( 1 - obj_optCART.criteria_weight ) .* ...
       ( ( 1 ./ N_tr + 1./ N_est ) * sum( var_combined( 1 , 1 , : ) ) );

end
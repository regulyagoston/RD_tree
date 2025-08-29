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

function [ sQ, e_tau_sq_1, e_tau_sq_2, avg_bias_sq, var_part ]  = treat_local_lbl_RDD_crit( leaf_nodes , obj_sample , obj_optCART )

[ tau_i, v_j, sigma2_bc_j , bias_j, h_j ] = collectInputs( leaf_nodes , true , ...
                                            'tau_i',...
                                            'v_j', 'var_tau_bc_j' ,...
                                            'bias_j', 'h_j' );

% Additional information to be used
N_tr  = obj_sample.n_tr;
N_est = obj_sample.n_est;
nPi   = numel( leaf_nodes );
%
n_j_tr = [leaf_nodes.n_j_tr]';
p_j_tr = n_j_tr/N_tr;
n_j_est = [leaf_nodes.n_j_est]';
p_j_est = n_j_est/N_est;

sigma2_bc_j = reshape( sigma2_bc_j, [nPi,1] );

% Expected value of SQUARED treatments
e_tau_sq_1 = mean( tau_i.^2 );
e_tau_sq_2 = sum( sigma2_bc_j .* p_j_tr );
e_tau_sq = -( e_tau_sq_1 - e_tau_sq_2 ); 
% Averaged bias square
avg_bias_sq = sum( h_j.^2 .* (bias_j.^2) .* p_j_est);
% Overall var scalred
sum_var = sum( v_j ./ h_j );

% Bandwidth these are the same below and above at the moment...
var_part = sum_var/N_est;



% Criterion
sQ = e_tau_sq + avg_bias_sq + var_part;

end
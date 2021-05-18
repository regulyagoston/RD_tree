%% Estimates the treatment effect using simple means
% 
% Input:
%   node - nodeProp outect
%       contains the node/leaf which needs to be estimated
%   data - sample outect
%       contains all the data for the estimation
%   optTree - optCART outect
%       contains relevant information about what to estimate
%   tr - 1x1 logical value
%       true  - it uses the training sample
%       false - it uses the estiation sample
%
% Output
%   out - structure
%       Y_hat_i   - predicted values
%       tau_i     - treamtent values 
%       sigma2_pj - variance of the treated units
%       sigma2_mj - variance of the untreated units

function out = est_treat_lbl_mean( node , data , optTree , tr )

%
out = struct;

% Honest and treatment type values
if optTree.honest && any( strcmp( optTree.type , {'CATE','CLATE'} ) )
    out.n_p_j_tr  = sum( node.logID_tr &  data.ITT_tr );
    out.n_m_j_tr  = sum( node.logID_tr & ~data.ITT_tr );
    out.n_p_j_est = sum( node.logID_est &  data.ITT_est );
    out.n_m_j_est = sum( node.logID_est & ~data.ITT_est );
    out.p_p_j_est = out.n_p_j_est ./ ( out.n_p_j_est + out.n_m_j_est );
    out.p_m_j_est = out.n_m_j_est ./ ( out.n_p_j_est + out.n_m_j_est );
end


%% Get the data to work with (training or estimation sample)
if tr
    [ Y , ~ , W , ITT ] = get_node_tr_Sample( node , data );
    n_p = out.n_p_j_tr;
    n_m = out.n_m_j_tr;
else
    [ Y , ~ , W , ITT ] = get_node_est_Sample( node , data );
    n_p = out.n_p_j_est;
    n_m = out.n_m_j_est;
end
n_ell = numel( Y );
nITT = ~ITT;

%% Leave-one-out estimator
if isprop( optTree , 'LOO' ) && optTree.LOO
    % Only calculates until one of the treated/untreated is done
    n_loo = min( n_p , n_m );
    % Slots
    d_p         = NaN( n_loo , 1 );
    d_m         = NaN( n_loo , 1 );
    tau_i       = NaN( n_loo , 1 );
    sigma2_pj_i = NaN( n_loo , 1 );
    sigma2_mj_i = NaN( n_loo , 1 );
    % Observations
    Yt = Y(  ITT );
    Yu = Y( nITT );
    for i = 1 : n_loo
        Yt_i = Yt( [ 1 : i - 1 , i + 1 : end ] , : );
        Yu_i = Yu( [ 1 : i - 1 , i + 1 : end ] , : );
        d_p( i ) = mean( Yt_i );
        d_m( i ) = mean( Yu_i );
        sigma2_pj_i( i ) = var( Yt_i );
        sigma2_mj_i( i ) = var( Yu_i );
        tau_i( i ) = d_p( 1 ) - d_m( 1 );
    end
    % Predicted values
    out.Y_hat_i       = NaN( n_ell , 1 );
    out.Y_hat_i(  W ) = mean( d_p );
    out.Y_hat_i( ~W ) = mean( d_m );
    % Treatment values
    out.tau_i     = repmat( mean( tau_i ) , [ n_ell , 1 ] );
    out.sigma2_p_j = mean( sigma2_pj_i );
    out.sigma2_m_j = mean( sigma2_mj_i );
else
    %% Simple sample-estimator
    d_p = mean( Y(  ITT ) );
    d_m = mean( Y( nITT ) );
    % Predicted values
    out.Y_hat_i       = NaN( n_ell , 1 );
    out.Y_hat_i(  ITT ) = d_p;
    out.Y_hat_i( nITT ) = d_m;
    % Treatment values
    out.tau_j = d_p - d_m;
    out.tau_i = repmat( out.tau_j , [ n_ell , 1 ] );
    out.sigma2_p_j = var( Y(  ITT ) );
    out.sigma2_m_j = var( Y( nITT ) );
    out.tau_se_j = sqrt( out.sigma2_p_j ./ n_p + out.sigma2_m_j ./ n_m );
end


end
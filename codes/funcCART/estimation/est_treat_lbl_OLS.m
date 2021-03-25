%% Estimates the treatment effect using simple OLS: y_+ = 1 + beta_1 * x_+ + beta_2 * x^2_+ + ... beta_p * x^p_+
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
%       tau_i      - treamtent values 
%       Mi_p_j     - scaled inverses of the covariate matrix, treated units homoskedastic case
%       Mi_m_j     - scaled inverses of the covariate matrix, untreated units homoskedastic case
%       sigma2_p_j - variance of the treated units
%       sigma2_m_j - variance of the untreated units
%       p_p_j      - share of treated units in the leaf
%       p_m_j      - share of untreated units in the leaf

function out = est_treat_lbl_OLS( node , data , optTree , tr )

%
out = struct;

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
end

fuzzy = strcmp( optTree.type , 'CLATE' ) && strcmp( optTree.criteria , 'RDD' );
%% Get the data to work with (training or estimation sample)
if tr
    [ Y , X , W , ITT , cl ] = get_node_tr_Sample( node , data );
    n_p = out.n_p_j_tr;
    n_m = out.n_m_j_tr;
else
    [ Y , X , W , ITT , cl ] = get_node_est_Sample( node , data );
    n_p = out.n_p_j_est;
    n_m = out.n_m_j_est;
end
% Number of observation in the leaf
n_ell = numel( Y );
% Not intended-to-treat observations
nITT = ~ITT;
% Create the regressors
p  = data.nCov;
X0 = [ ones( n_ell , 1 ) , X( : , 1 : p ) ];
if any( p + 1 >= n_p | p + 1 >= n_m )
    error('Not enough observation in the leaf! Impossible to calculate the treatment effect!')
end

%% Leave-one-out estimator
if isprop( optTree , 'LOO' ) && optTree.LOO
    error('Not updated...')
    % Only calculates until one of the treated/untreated is done
    n_loo = min( n_p , n_m );
    % Slots
    d_p         = NaN( n_loo , 1 );
    d_m         = NaN( n_loo , 1 );
    tau_i       = NaN( n_loo , 1 );
    sigma2_pj_i = NaN( n_loo , 1 );
    sigma2_mj_i = NaN( n_loo , 1 );
    Mi_p        = NaN( p + 1 , p + 1 , n_loo );
    Mi_m        = NaN( p + 1 , p + 1 , n_loo );
    % Observations
    Yt = Y(  ITT );
    Yu = Y( nITT );
    X0_p = X0(  ITT , : );
    X0_m = X0( nITT , : );
    for i = 1 : n_loo
        Yt_i = Yt( [ 1 : i - 1 , i + 1 : end ] , : );
        Yu_i = Yu( [ 1 : i - 1 , i + 1 : end ] , : );
        Xt_i = X0_p( [ 1 : i - 1 , i + 1 : end ] , : );
        Xu_i = X0_m( [ 1 : i - 1 , i + 1 : end ] , : );
        % treatment_-i
        d_p( i ) = lscov( Xt_i , Yt_i );
        d_m( i ) = lscov( Xu_i , Yu_i );
        tau_i( i ) = d_p( 1 ) - d_m( 1 );
        % sigma2_-i
        sigma2_pj_i( i ) = 1./( n_p - p - 2 ) .* sum( ( Yt_i - Xt_i * d_p( i ) ).^2 );
        sigma2_mj_i( i ) = 1./( n_m - p - 2 ) .* sum( ( Yu_i - Xu_i * d_m( i ) ).^2 );
        % Mi_-i
        Mi_p( : , : , i ) = inv( 1 ./ n_p .* ( Xt_i' * Xt_i ) );
        Mi_m( : , : , i ) = inv( 1 ./ n_m .* ( Xu_i' * Xu_i ) );
    end
    % Averages values
    out.tau_i      = repmat( mean( tau_i ) , [ n_ell , 1 ] );
    out.sigma2_p_j = mean( sigma2_pj_i );
    out.sigma2_m_j = mean( sigma2_mj_i );
    out.Mi_p_j = squeeze( mean( Mi_p , 3 ) );
    out.Mi_m_j = squeeze( mean( Mi_m , 3 ) );
else
    %% Simple OLS estimator for the treatment
    Y_p  = Y(  ITT );
    Y_m  = Y( nITT );
    X0_p = X0(  ITT , : );
    X0_m = X0( nITT , : );
    d_p = lscov( X0_p , Y_p );
    d_m = lscov( X0_m , Y_m );
    out.d_p = d_p;
    out.d_m = d_m;
    % Treatment values - diff in the intercept
    out.tau_i = repmat( d_p( 1 ) - d_m( 1 ) , [ n_ell , 1 ] );
    % Inverses for homoskedastic case
    out.Mi_p_j = inv( 1 ./ n_p .* ( X0_p' * X0_p ) );
    out.Mi_m_j = inv( 1 ./ n_m .* ( X0_m' * X0_m ) );
    
    if fuzzy
        W_p  = W(  ITT );
        W_m  = W( nITT );
        dw_p = lscov( X0_p , W_p );
        dw_m = lscov( X0_m , W_m );
        out.dw_p = dw_p;
        out.dw_m = dw_m;
        % First-stage Treatment values - diff in the intercept
        out.tau_y_i = out.tau_i;
        tau_w_i = dw_p( 1 ) - dw_m( 1 );
        out.tau_w_i = repmat( tau_w_i , [ n_ell , 1 ] );
        if out.tau_w_i < optTree.CLATE_crit_tau_w
            % This is going to make the EMSE Inf Large -> never choose it.
            out.sigma2_m_j = Inf( p + 1 , p + 1 );
            out.sigma2_p_j = Inf( p + 1 , p + 1 );
            out.tau_i      = zeros( n_ell , 1 );
            return;
        end
        % For variance estimation
        tau_y_i = out.tau_i( 1 );
        % Save actual LATE
        out.tau_i = out.tau_y_i ./ out.tau_w_i;
    else
        W_p = NaN;
        W_m = NaN;
        dw_p = NaN;
        dw_m = NaN;
        tau_w_i = NaN;
        tau_y_i = NaN;
    end
    
    % Different estimators for variance
    if ~isnan( cl( 1 ) )
        cl_p = cl(  ITT );
        cl_m = cl( nITT );
    else
        cl_p = NaN;
        cl_m = NaN;
    end
    
    if strcmp( optTree.criteria , 'RDD' )
        [ out.sigma2_m_j , out.var_y_m_j , out.var_w_m_j ] = get_var_estimator_RDD( ...
                                              Y_m , W_m , X0_m , d_m , dw_m , out.Mi_m_j , ...
                                              n_m , p , cl_m , ...
                                              fuzzy , tau_y_i , tau_w_i , optTree.varEstimator );
        [ out.sigma2_p_j , out.var_y_p_j , out.var_w_p_j ] = get_var_estimator_RDD( ...
                                              Y_p , W_p , X0_p , d_p , dw_p , out.Mi_p_j , ...
                                              n_p , p , cl_p , ...
                                              fuzzy , tau_y_i , tau_w_i , optTree.varEstimator );
        
        tau_y_var = out.var_y_m_j + out.var_y_p_j;
        out.tau_w_var = out.var_w_m_j + out.var_w_p_j;
    else
       error('Not implemented variance, only sharp RDD or fuzzy RDD') 
    end
    out.tau_var_j = tau_y_var( 1 , 1 );
    out.tau_se_j  = sqrt( out.tau_var_j );
end


end
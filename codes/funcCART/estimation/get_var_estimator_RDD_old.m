%% Different estimators used for variance in RDD
%
% Y_p  - outcome variable for above threshold (Nx1)
% W_p  - treatment variable for above threshold (Nx1)
% X0_p - covariates (constant+running variable's polinomial) for treated
%   (N x p+1)
% d_p   - estimated coefficients for X0_p ( p+1 x 1 ) on outcome
% dw_p  - estimated coefficients for X0_p ( p+1 x 1 ) on treatment
% Mi_p_j - inverse of ( X0_p' * X0_p ) for each leaf ( p+1 x p+1 x numLeaf )
%
% Same for control
%
% n_p - number of observations for treated units (scalar)
% n_c - number of observations for control units
% p   - number of covariates in X0_p
% varEstimator - type of variance estimator:
%       - 'simple' - homscedastic case
%       - 'hce-0', 'hce-1' - heteroscedastic robust se
%       - 'clustered' - clustered se
%
% cl_p - observations to be clustered for treated units
% cl_m - observations to be clustered for treated units



function [ Sigma_Y_m , Sigma_Y_p , ...
           V_Y_m , V_Y_p , var_tau_y , ...
           V_W_m , V_W_p , var_tau_w ] = get_var_estimator_RDD_old( Y_p , W_p , X0_p , d_p , dw_p , Mi_p_j , ...
                                                                Y_m , W_m , X0_m , d_m , dw_m , Mi_m_j , ...
                                                                n_p , n_m , p , cl_p , cl_m , ...
                                                                fuzzy , varEstimator)
%% Predict the outcomes for EQ(s)
% 1) jointly predict the outcomes and calculate the residuals
if fuzzy
    pred_m = X0_m * [ d_m , dw_m ];
    pred_p = X0_p * [ d_p , dw_p ];
    res_m = [ Y_m , W_m ] - pred_m;
    res_p = [ Y_p , W_p ] - pred_p;
else
    pred_m = X0_m * d_m;
    pred_p = X0_p * d_p;
    res_m  = Y_m - pred_m;
    res_p  = Y_p - pred_p;
end

% Adjust the inverse of X'X with the number of observations from the given
% (other) sample
Mi_p_j = 1 ./ n_p .* Mi_p_j;
Mi_m_j = 1 ./ n_m .* Mi_m_j;

% Add weights for different type of estimators
if any( strcmp( varEstimator , {'simple','hce-0'} ) )
    w_m1 = 1;
    w_p1 = 1;
elseif any( strcmp( varEstimator , {'hce-1','clustered'} ) )
    w_m1 = sqrt( n_m ./ ( n_m - p - 1 ) );
    w_p1 = sqrt( n_p ./ ( n_p - p - 1 ) );
else
    error('No such variance estimator is implemented for fuzzy RDD!')
end
res_m = w_m1 .* res_m;
res_p = w_p1 .* res_p;

% 3) Calculate the Var-Covar matrix:
if fuzzy
    % 3/a) Scaling vector for tau_y and tau_t
    tau_y = d_p( 1 )  - d_m( 1 );
    tau_w = dw_p( 1 ) - dw_m( 1 );
    s_Y = [ 1./ tau_w , -( tau_y./tau_w.^2 ) ];
    s_T = [ 0 , 1 ];
    d = 2;
else
    s_Y = 1;
    d = 1;
    % For required output
    V_W_m = NaN;
    V_W_p = NaN;
    var_tau_w = NaN;
end

% 3/b estimate sigma2 for all regressions
% number of covariates
k = size( X0_m , 2 );
Sigma_Y_m = zeros( k , k );
Sigma_Y_p = zeros( k , k );
if fuzzy
    Sigma_W_m = zeros( k , k );
    Sigma_W_p = zeros( k , k );
end
ng_log = nargout > 2 && fuzzy;

if strcmp( varEstimator , 'clustered' )
    clu_m = unique( cl_m , 'stable' );
    clu_p = unique( cl_p , 'stable' );
    g_m = numel( clu_m );
    g_p = numel( clu_p );
    w_m = ( n_m - 1 ) ./ ( n_m - p - 1  ) .* ( g_m ./ ( g_m - 1 ) );
    w_p = ( n_p - 1 ) ./ ( n_p - p - 1  ) .* ( g_p ./ ( g_p - 1 ) );
    
    for i = 1 : max( g_m , g_p )
        if i <= g_m
            ind_m = cl_m == clu_m( i );
            Xi = X0_m( ind_m , : );
            ri = res_m( ind_m, : );
            for l = 1 : d
                for j = 1 : d
                    Sigma_Y_m = Sigma_Y_m + ( ( Xi' * ( s_Y( l ) .* ri( : , l ) ) ) * ( Xi' * ( s_Y( j ) .* ri( : , j ) ) )' )';
                    if ng_log
                        Sigma_W_m = Sigma_W_m + ( ( Xi' * ( s_T( l ) .* ri( : , l ) ) ) * ( Xi' * ( s_T( j ) .* ri( : , j ) ) )' )';
                    end
                end
            end
        end
        if i <= g_p
            ind_p = cl_p == clu_p( i );
            Xi = X0_p( ind_p , : );
            ri = res_p( ind_p, : );
            for l = 1 : d
                for j = 1 : d
                    Sigma_Y_p = Sigma_Y_p + ( ( Xi' * ( s_Y( l ) .* ri( : , l ) ) ) * ( Xi' * ( s_Y( j ) .* ri( : , j ) ) )' )';
                    if ng_log
                        Sigma_W_p = Sigma_W_p + ( ( Xi' * ( s_T( l ) .* ri( : , l ) ) ) * ( Xi' * ( s_T( j ) .* ri( : , j ) ) )' )';
                    end
                end
            end
        end
    end
    % Scale it with w for later usage in the criterion
    Sigma_Y_p = Sigma_Y_p .* w_p;
    Sigma_Y_m = Sigma_Y_m .* w_m;
    if ng_log
        Sigma_W_p = Sigma_W_p .* w_p;
        Sigma_W_m = Sigma_W_m .* w_m;
    end
elseif any( strcmp( varEstimator , {'hce-0','hce-1'} ) )
    for i = 1 : d
        % Squared term then product of the two parts
        SS_m = res_m( : , i ) .* res_m;
        SS_p = res_p( : , i ) .* res_p;
        for j = 1 : d
            % Calculate individually Imbens-Lemieux p. 630
            Sigma_Y_m = Sigma_Y_m + ( X0_m .* ( s_Y( i ) .* s_Y( j ) ) .* SS_m( : , j ) )' * X0_m;
            Sigma_Y_p = Sigma_Y_p + ( X0_p .* ( s_Y( i ) .* s_Y( j ) ) .* SS_p( : , j ) )' * X0_p;
            if ng_log
                Sigma_W_m = Sigma_W_m + ( X0_m .* ( s_T( i ) .* s_T( j ) ) .* SS_m( : , j ) )' * X0_m;
                Sigma_W_p = Sigma_W_p + ( X0_p .* ( s_T( i ) .* s_T( j ) ) .* SS_p( : , j ) )' * X0_p;
            end
        end
    end 
else
    % Squared term then product of the two parts
    SS_m = 1 ./ ( n_m - p - 1 ) .* res_m' * res_m;
    SS_p = 1 ./ ( n_p - p - 1 ) .* res_p' * res_p;
    Sigma_Y_m = 0;
    Sigma_Y_p = 0;
    Sigma_W_m = 0;
    Sigma_W_p = 0;
    for i = 1 : d
        for j = 1 : d
            % Calculate individually Imbens-Lemieux p. 630
            Sigma_Y_m = Sigma_Y_m + ( s_Y( i ) .* s_Y( j ) ) .* SS_m( i , j );
            Sigma_Y_p = Sigma_Y_p + ( s_Y( i ) .* s_Y( j ) ) .* SS_p( i , j );
            if ng_log
                Sigma_W_m = Sigma_W_m + ( s_T( i ) .* s_T( j ) ) .* SS_m( i , j );
                Sigma_W_p = Sigma_W_p + ( s_T( i ) .* s_T( j ) ) .* SS_p( i , j );
            end
        end
    end 
end

if ~strcmp( varEstimator , 'simple' )
    % Scale it for EMSE criterion
    Sigma_Y_m = Sigma_Y_m ./ n_m;
    Sigma_Y_p = Sigma_Y_p ./ n_p;
end

if nargout > 2
    if any( strcmp( varEstimator , {'hce-0','hce-1','clustered'} ) )
        % Variance below and above in sandwich form with the scaling
        V_Y_m = Mi_m_j * ( n_m .* Sigma_Y_m ) * Mi_m_j;
        V_Y_p = Mi_p_j * ( n_p .* Sigma_Y_p ) * Mi_p_j;
        if ng_log
            V_W_m = Mi_m_j * Sigma_W_m * Mi_m_j;
            V_W_p = Mi_p_j * Sigma_W_p * Mi_p_j;
            Sigma_W_m = Sigma_W_m ./ n_m;
            Sigma_W_p = Sigma_W_p ./ n_p;
        end
    else
        V_Y_m = Sigma_Y_m .* Mi_m_j;
        V_Y_p = Sigma_Y_p .* Mi_p_j;
        if ng_log
            V_W_m = Sigma_W_m .* Mi_m_j;
            V_W_p = Sigma_W_p .* Mi_p_j;
        end
    end

if nargout > 4
    % The overall var-covar matrix
    var_tau_y = V_Y_m + V_Y_p;
    var_tau_w = V_W_m + V_W_p;
end

end
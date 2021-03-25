%% Different estimators used for variance in RDD
%
% Y  - outcome variable (Nx1)
% W  - treatment variable (Nx1) - can be NaN
% X0 - covariates (constant+running variable's polinomial) (N x p+1)
% d_y - estimated coefficients for X0 ( p+1 x 1 ) on outcome
% d_w - estimated coefficients for X0 ( p+1 x 1 ) on treatment - can be NaN
% Mi  - inverse of ( X0' * X0 ) for each leaf ( p+1 x p+1 x numLeaf )
%
% n - number of observations
% p - number of covariates in X0_p
% varEstimator - type of variance estimator:
%       - 'simple' - homscedastic case
%       - 'hce-0', 'hce-1' - heteroscedastic robust se
% cl - clusters for observations



function [ Sigma_Y , V_Y , V_W ] = get_var_estimator_RDD( Y , W , X0 , d_y , d_w , Mi , ...
                                                          n , p , cl , fuzzy , tau_y , tau_w , varEstimator)
%% Predict the outcomes for EQ(s)
% 1) jointly predict the outcomes and calculate the residuals
if fuzzy
    pred = X0 * [ d_y , d_w ];
    res = [ Y , W ] - pred;
else
    pred = X0 * d_y;
    res  = Y - pred;
end

% Adjust the inverse of X'X with the number of observations from the given
% (other) sample
Mi = 1 ./ n .* Mi;

% Select the proper variance estimator
ls = false;
l0 = false;
l1 = false;
switch varEstimator
    case 'simple'
        ls = true;
    case 'hce-0'
        l0 = true;
    case 'hce-1'
        l1 = true;
end

% Decide if there are clusters
l_cl = all( ~isnan( cl ) );

% Decide if there is need to compute the variance of the treatment eq
l_teq = nargout > 2 && fuzzy;

% Add weights for different type of estimators
if ls || l0
    w1 = 1;
elseif l1
    w1 = sqrt( n ./ ( n - p - 1 ) );
else
    error('No such variance estimator is implemented for fuzzy RDD!')
end
res = w1 .* res;

% 3) Calculate the Var-Covar matrix:
if fuzzy
    % 3/a) Scaling vector for tau_y and tau_t
    s_Y = [ 1./ tau_w , -( tau_y./tau_w.^2 ) ];
    s_T = [ 0 , 1 ];
    d = 2;
else
    s_Y = 1;
    d = 1;
    % For required output
    V_W = NaN;
end

% 3/b estimate sigma2 for all regressions
% number of covariates
k = size( X0 , 2 );
Sigma_Y = zeros( k , k );
if fuzzy
    Sigma_W = zeros( k , k );
end


if l_cl
    % Clustered SE's
    clu = unique( cl );
    g = numel( clu );
    w2 = ( n - 1 ) ./ ( n - p - 1  ) .* ( g ./ ( g - 1 ) );
    for i = 1 : g
        ind = cl == clu( i );
        Xi = X0( ind , : );
        ri = res( ind, : );
        for l = 1 : d
            for j = 1 : d
                Sigma_Y = Sigma_Y + ( ( Xi' * ( s_Y( l ) .* ri( : , l ) ) )...
                                    * ( Xi' * ( s_Y( j ) .* ri( : , j ) ) )' )';
                if l_teq
                    Sigma_W = Sigma_W + ( ( Xi' * ( s_T( l ) .* ri( : , l ) ) ) ...
                                        * ( Xi' * ( s_T( j ) .* ri( : , j ) ) )' )';
                end
            end
        end
    end
    % Scale it with w for later usage in the criterion
    Sigma_Y = Sigma_Y .* w2;
    if l_teq
        Sigma_W = Sigma_W .* w2;
    end
elseif l0 || l1
    % HCE-0 and HCE-1 estimators
    for i = 1 : d
        % Squared term then product of the two parts
        SS = res( : , i ) .* res;
        for j = 1 : d
            % Calculate individually Imbens-Lemieux p. 630
            Sigma_Y = Sigma_Y + ( X0 .* ( s_Y( i ) .* s_Y( j ) ) .* SS( : , j ) )' * X0;
            if l_teq
                Sigma_W = Sigma_W + ( X0 .* ( s_T( i ) .* s_T( j ) ) .* SS( : , j ) )' * X0;
            end
        end
    end 
else
    % Homoscedastic case
    % Squared term then product of the two parts
    SS = 1 ./ ( n - p - 1 ) .* res' * res;
    Sigma_Y = 0;
    Sigma_W = 0;
    for i = 1 : d
        for j = 1 : d
            % Calculate individually Imbens-Lemieux p. 630
            Sigma_Y = Sigma_Y + ( s_Y( i ) .* s_Y( j ) ) .* SS( i , j );
            if l_teq
                Sigma_W = Sigma_W + ( s_T( i ) .* s_T( j ) ) .* SS( i , j );
            end
        end
    end 
end

if ~ls
    % Scale it for EMSE criterion
    Sigma_Y = Sigma_Y ./ n;
end

if nargout > 2
    if l0 || l1
        % Variance below and above in sandwich form with the scaling
        V_Y = Mi * ( Sigma_Y .* n ) * Mi;
        if l_teq
            V_W = Mi * Sigma_W * Mi;
            %Sigma_W = Sigma_W ./ n;
        end
    else
        V_Y = Sigma_Y .* Mi;
        if l_teq
            V_W = Sigma_W .* Mi;
        end
    end
end

end
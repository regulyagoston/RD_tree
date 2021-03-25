%% CT-Athey-Imbens splitting criterion - implementation of C++ code
%
%
% to-do: there is a slight difference between updating formula and simple
% OLS. I do not know why, but sometimes there is a 10e-3 difference, which
% can be considered as substantial!

function split  = split_CT_RDD( Y , W , Z , X , ITT , cl , Z_est , ITT_est , optTree )


%% Initialization
n       = numel( Y );
n_t     = sum( ITT );
n_c     = n - n_t;
p       = size( X , 2 );
alpha   = optTree.criteria_weight;
minObs  = optTree.minObs;
% Whether it is fuzzy design or not
fuzzy = strcmp( optTree.type , 'CLATE' );

%% First calculate the node's variances
X0   = [ ones( n , 1 ) , X ];
uITT = ~ITT;
uITT_est    = ~ITT_est;
Y_t  = Y(  ITT );
Y_c  = Y( uITT );
X0_t = X0(  ITT , : );
X0_c = X0( uITT , : );
cl_t = cl(  ITT );
cl_c = cl( uITT );

if fuzzy
    W_t = W(  ITT );
    W_c = W( uITT );
else
    W_t = NaN( n_t , 1 );
    W_c = NaN( n_c , 1 );
end


%%%%%%%%%%
%% Check for each splitting candidate
%%%%%%%%%%
% To be efficient and fast only check those which satisfies:
%   there are enough observation in each split (lower and upper parts)
%   there are enough observation in each bucket (number of treated and control)
%
% I use Sherman-Morrison formula with recursive updating:
%   faster than lscov even if it needs to compute for each observation

%% First check for enough observations in the 'tails'
left_tr     = cumsum(  ITT );
left_con    = cumsum( uITT );
left_tr_e   = cumsum(  ITT_est );
left_con_e  = cumsum( uITT_est );
right_tr    = n_t - left_tr;
right_con   = n_c - left_con;
right_tr_e  = sum(  ITT_est ) - left_tr_e;
right_con_e = sum( uITT_est ) - left_con_e;
% Logicals to have enough observations on both side in both samples
logNumObs   = left_tr >= minObs & ...                         
              left_con >= minObs & ...                        
              right_tr >= minObs & ...                        
              right_con >= minObs;
logNumObs_e = left_tr_e >= minObs & ...                       
              left_con_e >= minObs & ...                      
              right_tr_e >= minObs & ...               
              right_con_e >= minObs;
Z_e_valid = Z_est( logNumObs_e );

if ~isempty( Z_e_valid )
    logNumObs = logNumObs & ...
            Z >= Z_e_valid( 1 ) & ...
            Z <= Z_e_valid( end );
else
    split = NaN;
    return;
end
               
if ~any( logNumObs & ITT  ) || ~any( logNumObs & uITT  )
    split = NaN;
    return;
end

% Get the first and last valid observations to iterate through
aux       = Z( logNumObs );
lNumObsZ  = Z >= aux( 1 ) & Z <= aux( end );
v_obs_tr  = Z( lNumObsZ &  ITT );
v_obs_co  = Z( lNumObsZ & uITT );
% Intent-to-treat values and controls
Ztr = Z( ITT , 1 );
Zco = Z( uITT , 1 );
% where does valid treated and control observations would start
f_obs_tr_aux  = max( find( v_obs_tr( 1 )   == Ztr , 1 , 'last' ) , minObs );
l_obs_tr_aux  = min( find( v_obs_tr( end ) == Ztr , 1 , 'first' ) , n - minObs );
f_obs_co_aux  = max( find( v_obs_co( 1 )   == Zco , 1 , 'last' ) , minObs );
l_obs_co_aux  = min( find( v_obs_co( end ) == Zco , 1 , 'first' ) , n - minObs );
% Get the indexes for overall observations
f_obs_all = find( Z == max( Ztr( f_obs_tr_aux ) , Zco( f_obs_co_aux ) ) , 1 , 'last' );
l_obs_all = find( Z == min( Ztr( l_obs_tr_aux ) , Zco( l_obs_co_aux ) ) , 1 , 'first' );
% Adjust the first and last observations for both tr and co
f_obs_tr = find( Ztr <= Z( f_obs_all ) , 1 , 'last' );
f_obs_co = find( Zco <= Z( f_obs_all ) , 1 , 'last' );
l_obs_tr = find( Ztr >= Z( l_obs_all ) , 1 , 'first' );
l_obs_co = find( Zco >= Z( l_obs_all ) , 1 , 'first' );

if any( f_obs_all >= l_obs_all || f_obs_tr >= l_obs_tr || f_obs_co >= l_obs_co )
    split = NaN;
    return;
end

% Let compute the parameter values for all the possible split values
n_v = l_obs_all - f_obs_all + 1;


% But we have to deal with the same values in Z:
% Number of iterations: unique values in valid Z
% All observations
[ Zuq , idZ ] = unique( Z( f_obs_all : l_obs_all ) );
% The starting indeces of the valid Z values:
idZ_all_uc = idZ + f_obs_all - 2;
% The first valid value must be dropped, while the algorithm only considers
%   values which are less or equal to a certain value
idZ_all = idZ_all_uc( 2 : end );
% Number of valid splits
n_v_sc = numel( idZ_all );

% Valid splits among ITT and control group
l_aux_tr = Ztr >= Zuq( 1 );
f_obs_tr_v = find( l_aux_tr , 1 , 'first' );
[ Zuq_tr , idZ_tr ] = unique( Ztr( l_aux_tr & Ztr <= Zuq( end ) ) );
l_aux_co = Zco >= Zuq( 1 );
f_obs_co_v = find( l_aux_co , 1 , 'first' );
[ Zuq_co , idZ_co ] = unique( Zco( l_aux_co & Zco <= Zuq( end ) ) );

% Valid and unique values' indexes for treated and control units based on
% each variable
idZ_tr_var = idZ_tr( 2 : end ) - 1;
idZ_co_var = idZ_co( 2 : end ) - 1;

% If there are no unique treated/control variables, return
if any( isempty( idZ_tr_var ) || isempty( idZ_co_var ) )
    split = NaN;
    return;
end

% %%
% Short-cut:
% If there are many same values in Z, one can speed up the process with
% simple OLS instead of updating algorithm, where much less computation is needed.
%   experience: if there are less than 100,000 observations 3 times more
%   estimation is approx the same computation time. Above that - because of
%   the inversion - 10 times is better approx (of course it is an exponential
%   fnc between computation time and number of obs...)

if ( n_v_sc * 3 < n_v && n_v < 100000 ) || ( n_v_sc * 10 < n_v )
    
    % Simple least squares:
    %   Estimate poi-s for each unique values
    [ b_y_t_a , b_y_t_b , b_y_c_a , b_y_c_b , ...
      Mi_t_a , Mi_t_b , Mi_c_a , Mi_c_b , ...
      b_w_t_a , b_w_t_b , b_w_c_a , b_w_c_b ] = ...
                               ls_est_coeffs( p , n_t , n_c , idZ_tr_var , idZ_co_var , ...
                                              X0_t , X0_c , Y_t , Y_c , W_t , W_c , fuzzy );
    
else 
    % Use the Sherman-Morrison formula with updating the values
    
    % 1) Get the number of updates for treated and non-treated units
    n_v_t     = l_obs_tr - f_obs_tr + 1;
    n_v_c     = l_obs_co - f_obs_co + 1;
    
    % 2) estimate all the poi-s regardless of replication
    [ b_y_t_a_r , b_y_t_b_r , b_y_c_a_r , b_y_c_b_r ,...
      Mi_t_a_r , Mi_t_b_r , Mi_c_a_r , Mi_c_b_r , ...
      b_w_t_a_r , b_w_t_b_r , b_w_c_a_r , b_w_c_b_r ] = ...
               sm_est_coeffs( p , n_t , n_c , n_v_t , n_v_c , f_obs_tr , l_obs_tr , f_obs_co , l_obs_co , ...
                              X0_t , X0_c , Y_t , Y_c , W_t , W_c , fuzzy );
    
    % 3) remove poi-s, where Ztr or Zco does not change:
    %  a) find the indexes
    [ ~ , idZ_tr_u ] = unique( Ztr );
    idZ_tr_u = idZ_tr_u( idZ_tr_u > f_obs_tr & idZ_tr_u <= l_obs_tr ) - f_obs_tr;
    [ ~ , idZ_co_u ] = unique( Zco );
    idZ_co_u = idZ_co_u( idZ_co_u > f_obs_co & idZ_co_u <= l_obs_co ) - f_obs_co;
    
    %   b) retain only those values:
    % for t above idZ_tr+1 because of Less or Equal rule
    % for t below idZ_tr
    b_y_t_a = b_y_t_a_r( idZ_tr_u + 1 , : );
    b_y_t_b = b_y_t_b_r( idZ_tr_u     , : );
    b_y_c_a = b_y_c_a_r( idZ_co_u + 1 , : );
    b_y_c_b = b_y_c_b_r( idZ_co_u     , : );
    Mi_t_a  = Mi_t_a_r( : , : , idZ_tr_u + 1 );
    Mi_t_b  = Mi_t_b_r( : , : , idZ_tr_u     );
    Mi_c_a  = Mi_c_a_r( : , : , idZ_co_u + 1 );
    Mi_c_b = Mi_c_b_r( : , : , idZ_co_u     );
    b_w_t_a = b_w_t_a_r( idZ_tr_u + 1 , : );
    b_w_t_b = b_w_t_b_r( idZ_tr_u     , : );
    b_w_c_a = b_w_c_a_r( idZ_co_u + 1 , : );
    b_w_c_b = b_w_c_b_r( idZ_co_u     , : );
    
end


%% Calculate honest splitting criterion
split_crit = NaN( n_v_sc , 1 );
n_v_t = max( numel( Zuq_tr ) - 1 , 1 );
n_v_c = max( numel( Zuq_co ) - 1 , 1 );
% Number of units for each valid splits
n_sc    = [ diff( idZ_all_uc )  ; n - idZ_all_uc( end ) ];
n_tr_sc = [ diff( idZ_tr ) ; n_t - idZ_tr( end ) + 1 ];
n_co_sc = [ diff( idZ_co ) ; n_c - idZ_co( end ) + 1 ];
n_sc( 1 ) = n_sc( 1 ) + f_obs_all - 1;
n_tr_sc( 1 ) = n_tr_sc( 1 ) + f_obs_tr - 1;
n_co_sc( 1 ) = n_co_sc( 1 ) + f_obs_co - 1;
% Because it is a split criteria, take the cumsum
n_sc_b    = cumsum( n_sc );
n_tr_sc_b = cumsum( n_tr_sc );

% Iteration number for treated and control for particular split:
%   Zuq is set such that it is valid for both treated and control group
ct_t = split_crit;
ct_c = split_crit;
ct_t( 1 ) = 1;
ct_c( 1 ) = 1;
t_count = true;
c_count = true;

% Only for checking the algorithm
tau_w_b = split_crit; tau_w_a = split_crit;
tau_b = split_crit; tau_a = split_crit;
var_b = NaN( p + 1 , p + 1 , numel( split_crit ) ); 
var_a = var_b;
logVar = any( strcmp( optTree.varEstimator , {'hce-0','hce-1','clustered'} ) );
s2Y_t_b = NaN;
s2Y_c_b = NaN;
for i = 1 : n_v_sc
    
    % Calculate the CATE/CLATE
    tau_y_b = b_y_t_b( ct_t( i ) , 1 ) - b_y_c_b( ct_c( i ) , 1 );
    tau_y_a = b_y_t_a( ct_t( i ) , 1 ) - b_y_c_a( ct_c( i ) , 1 );
    if fuzzy
        tau_w_b(i) = b_w_t_b( ct_t( i ) , 1 ) - b_w_c_b( ct_c( i ) , 1 );
        tau_w_a(i) = b_w_t_a( ct_t( i ) , 1 ) - b_w_c_a( ct_c( i ) , 1 );
        
        if ( tau_w_b( i ) > optTree.CLATE_crit_tau_w ) && ( tau_w_a( i ) > optTree.CLATE_crit_tau_w )
            tau_b(i) = tau_y_b ./ tau_w_b( i );
            tau_a(i) = tau_y_a ./ tau_w_a( i );
        else
            tau_b(i) = NaN;
            tau_a(i) = NaN;
            % In case it violates the CLATE restriction we need to
            % recompute the sigma2 in the next valid iteration!
            s2Y_t_b = NaN;
            s2Y_c_b = NaN;
        end
    else
        tau_b(i) = tau_y_b;
        tau_a(i) = tau_y_a;
    end
    if ~any( isnan( [ tau_b(i) , tau_a(i) ] ) )
        % Calculate the share of treated and untreated units
        n_b   = n_sc_b( i );
        n_t_b = n_tr_sc_b( ct_t( i ) );
        p_t_b = n_t_b ./ n_b;
        p_c_b = 1 - p_t_b;
        n_a   = n - n_b;
        n_t_a = n_t - n_t_b;
        p_t_a = n_t_a ./ n_a;
        p_c_a = 1 - p_t_a;
        % Calculate the variances:
        %   Y_p treated - below
        %   Y_m control - below
        
        % Current index for the variables
        idx_t = idZ_tr( ct_t( i ) ) + f_obs_tr_v - 1;
        idx_c = idZ_co( ct_c( i ) ) + f_obs_co_v - 1;
        
        % Sigma2 for treated units if there are new observations
        if isnan( s2Y_t_b( 1 ) ) || ct_t( i ) ~= ct_t( i - 1 )        
            s2Y_t_b = get_var_estimator_RDD( Y_t( 1 : idx_t , : ) , W_t( 1 : idx_t , : ) , X0_t( 1 : idx_t , : ) , ...
                                                b_y_t_b( ct_t( i ) , : )' , b_w_t_b( ct_t( i ) , : )' , Mi_t_b( : , : , ct_t( i ) ) , ...
                                                n_t_b , p , cl_t( 1 : idx_t , : ) , fuzzy , tau_y_b , tau_w_b( i ) , optTree.varEstimator );
            s2Y_t_a = get_var_estimator_RDD( Y_t( idx_t + 1 : end , : ) , W_t( idx_t + 1 : end , : ) , X0_t( idx_t + 1 : end , : ) , ...
                                                b_y_t_a( ct_t( i ) , : )' , b_w_t_a( ct_t( i ) , : )' , Mi_t_a( : , : , ct_t( i ) ) , ...
                                                n_t_a , p , cl_t( idx_t + 1 : end , : ) , fuzzy , tau_y_a , tau_w_a( i ) , optTree.varEstimator );
        end 
        if isnan( s2Y_c_b( 1 ) ) || ct_c( i ) ~= ct_c( i - 1 )        
            s2Y_c_b = get_var_estimator_RDD( Y_c( 1 : idx_c , : ) , W_c( 1 : idx_c , : ) , X0_c( 1 : idx_c , : ) , ...
                                                b_y_c_b( ct_c( i ) , : )' , b_w_c_b( ct_c( i ) , : )' , Mi_c_b( : , : , ct_c( i ) ) , ...
                                                n_b - n_t_b , p , cl_c( 1 : idx_c , : ) , fuzzy , tau_y_b , tau_w_b( i ) , optTree.varEstimator );
            s2Y_c_a = get_var_estimator_RDD( Y_c( idx_c + 1 : end , : ) , W_c( idx_c + 1 : end , : ) , X0_c( idx_c + 1 : end , : ) , ...
                                                b_y_c_a( ct_c( i ) , : )' , b_w_c_a( ct_c( i ) , : )' , Mi_c_a( : , : , ct_c( i ) ) , ...
                                                n_a - n_t_a , p , cl_c( idx_c + 1 : end , : ) , fuzzy , tau_y_a , tau_w_a( i ) , optTree.varEstimator );
        end
        
        % Decide if it is not heteroscedastic or yes
        if logVar
            var_b( : , : , i ) = ( Mi_t_b( : , : , ct_t( i ) ) * s2Y_t_b * Mi_t_b( : , : , ct_t( i ) ) ) ./ p_t_b ...
                                +( Mi_c_b( : , : , ct_c( i ) ) * s2Y_c_b * Mi_c_b( : , : , ct_c( i ) ) ) ./ p_c_b;
            var_a( : , : , i ) = ( Mi_t_a( : , : , ct_t( i ) ) * s2Y_t_a * Mi_t_a( : , : , ct_t( i ) ) ) ./ p_t_a ...
                                +( Mi_c_a( : , : , ct_c( i ) ) * s2Y_c_a * Mi_c_a( : , : , ct_c( i ) ) ) ./ p_c_a;
        else
            var_b( : , : , i ) = ( s2Y_t_b .* Mi_t_b( : , : , ct_t( i ) ) ) ./ p_t_b ...
                                +( s2Y_c_b .* Mi_c_b( : , : , ct_c( i ) ) ) ./ p_c_b;
            var_a( : , : , i ) = ( s2Y_t_a .* Mi_t_a( : , : , ct_t( i ) ) ) ./ p_t_a ...
                                +( s2Y_c_a .* Mi_c_a( : , : , ct_c( i ) ) ) ./ p_c_a;
        end
        
        % Split criterion
        split_crit( i ) = - alpha .* ( tau_b(i) .^ 2 .* n_b + tau_a(i) .^ 2 .* n_a ) ...
            + ( 1 - alpha ) .* ( squeeze( var_b( 1 , 1 , i ) + var_a( 1 , 1 , i ) ) );
    end

    % If there are new parameters
    if t_count && Zuq( i + 1 ) > Zuq_tr( ct_t( i ) )
        ct_t( i + 1 ) = ct_t( i ) + 1;
        if ct_t( i + 1 ) >= n_v_t
            ct_t( i + 1 ) = n_v_t;
            t_count = false;
        else
            while Zuq( i + 1 ) > Zuq_tr( ct_t( i + 1 ) )
                ct_t( i + 1 ) = ct_t( i + 1 ) + 1;
                if ct_t( i + 1 ) >= n_v_t
                    ct_t( i + 1 ) = n_v_t;
                    t_count = false;
                    break;
                end
            end
        end
    else
        ct_t( i + 1 ) = ct_t( i );
    end
    if c_count && Zuq( i + 1 ) > Zuq_co( ct_c( i ) )
        ct_c( i + 1 ) = ct_c( i ) + 1;
        if ct_c( i + 1 ) >= n_v_c
            ct_c( i + 1 ) = n_v_c;
            c_count = false;
        else
            while Zuq( i + 1 ) > Zuq_co( ct_c( i + 1 ) )
                ct_c( i + 1 ) = ct_c( i + 1 ) + 1;
                if ct_c( i + 1 ) >= n_v_c
                    ct_c( i + 1 ) = n_v_c;
                    c_count = false;
                    break;
                end
            end
        end
    else
        ct_c( i + 1 ) = ct_c( i );
    end
end

%% Use buckets for more robust criterion selection:
%   Athey-Imbens (2016)
% Check for possible NaN values
logValid = ~isnan( split_crit );

% If there are any candidates
if any( logValid )

    % Now sort them into buckets to smooth irregularities by individual
    % observations -> use buckets
    % Note Athey et al (2016) uses the buckets beforehand. However it does
    % not matter to do it in this way: with recursive updating I need to
    % calculate everything anyway...
    %
    % Find the proper number of observations in each buckets
    n_v_obs = sum( Z >= Zuq( 1 ) & Z <= Zuq( end ) );
    n_v_tr  = sum( Ztr >= Zuq( 1 ) & Ztr <= Zuq( end ) );
    
    % use for treated and untreated as well
    bucketTmp      = min( round( n_v_tr ./ optTree.obsBucket ) , ...
                          round( ( n_v_obs - n_v_tr ) ./ optTree.obsBucket ) ) .* 2;
    num_buckets    = min( bucketTmp , optTree.maxBucket );%max( minObs , min( bucketTmp , optTree.maxBucket ) );
    obs_in_buckets = round( n_v_obs ./ num_buckets );
    % Runs through the number of observations and select the end of each
    % bucket in the covariate
    logFin = false( n_v_sc , 1 );
    buc_num_i = NaN( n_v_sc , 1 );
    ct = 1;
    ct_t = 1;
    ct_c = 1;
    n_tr  = 0;
    n_con = 0;
    for i = 1 : n_v_sc
        % Count the number of treatment or control observation
        if Zuq( i ) >= Zuq_tr( ct_t )
            n_tr = n_tr + n_tr_sc( ct_t );
            if ct_t < n_v_t
                ct_t = ct_t + 1;
            end
        end
        if Zuq( i ) >= Zuq_co( ct_c )
            n_con = n_con + n_co_sc( ct_c );
            if ct_c < n_v_c
                ct_c = ct_c + 1;
            end
        end
        buc_num_i( i ) = ct;
        % If we have enough observations for both treated and control
        % than we accept as a candidate for the splitting value and restart
        % counting
        if n_tr >= obs_in_buckets && n_con >= obs_in_buckets
            logFin( i ) = true;
            n_tr = 0;
            n_con = 0;
            ct = ct + 1;
        end
    end
    % Finally select splitting criterion which satisfies the bucketing
    if any( logFin )
        % Find the best ID
        [ min_crit , mID ] = min( split_crit( logFin ) );
        % Find the largest values in the buckets
        % Exclude the last value while LEQ
        Zuq_a = Zuq( 1 : end - 1 );
        lZ_buckets = Zuq_a( logFin );
        % Splitting value is where the criterion is the smallest
        split = lZ_buckets( mID );
        % Alternatively Athey et al takes the averages - but this is stupid
        % Get the last value from the bucket for treatment and control
        %             b_tr  = aux_1( b_bucket == buc_num_i &  can_tr );
        %             b_con = aux_1( b_bucket == buc_num_i & ~can_tr );
        %             split = ( b_tr( end ) + b_con( end ) ) ./ 2;
    else
        split = NaN;
    end
else
    split = NaN;
end

end
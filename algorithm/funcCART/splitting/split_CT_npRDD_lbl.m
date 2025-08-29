%% CT- nonparametric splitting criterion with leaf-by-leaf bandwidths
%
%

function split  = split_CT_npRDD_lbl( Y , W , Z , X , ITT , cl , Z_est , ITT_est , X_est, c, N_tr, N_est, optTree )


%% Initialization
minObs  = optTree.minObs;
% Whether it is fuzzy design or not
fuzzy = strcmp( optTree.type , 'CLATE' );
uITT = ~ITT;
uITT_est    = ~ITT_est;
n_j = numel( Y );
n_t = sum( ITT );
n_c = sum( uITT );
Y_t  = Y(  ITT );
Y_c  = Y( uITT );
X_t = X(  ITT , : );
X_c = X( uITT , : );
cl_t = cl(  ITT );
cl_c = cl( uITT );

if fuzzy
    W_t = W(  ITT );
    W_c = W( uITT );
else
    W = NaN( n_j, 1);
end


%%%%%%%%%%
%% Check for each splitting candidate
%%%%%%%%%%
% To be efficient and fast only check those which satisfies:
%   there are enough observation in each split (lower and upper parts)
%   there are enough observation in each bucket (number of treated and control)
%

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
l_obs_tr_aux  = min( find( v_obs_tr( end ) == Ztr , 1 , 'first' ) , n_j - minObs );
f_obs_co_aux  = max( find( v_obs_co( 1 )   == Zco , 1 , 'last' ) , minObs );
l_obs_co_aux  = min( find( v_obs_co( end ) == Zco , 1 , 'first' ) , n_j - minObs );
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

if ( numel( Zuq ) > optTree.numSplit )
    % Get frequencies
    freq_Z = histcounts( Z( f_obs_all : l_obs_all ), Zuq );
    cum_prob = cumsum( freq_Z ) ./ numel( Z( f_obs_all : l_obs_all ) );
    % Requested percentages
    th_pctg = cumsum( ones( optTree.numSplit, 1 ) ) ./ optTree.numSplit;
    Zuq_iter = NaN( optTree.numSplit, 1 );
    for i = 1:optTree.numSplit
        [~,idx] = min( abs( cum_prob - th_pctg( i ) ) );
        Zuq_iter( i ) = Zuq( idx );
    end
    % In case of closest pctgs remains same
    Zuq_iter = unique( Zuq_iter );
else 
    Zuq_iter = Zuq;
end



%% Calculate the splitting criterion
n_split = numel( Zuq_iter );
split_crit = nan(n_split,1);
exit_flag = repmat({'success'},[n_split,1]);
kernFnc = str2func( [ optTree.kernel , '2' ] );

for i = 1:n_split
    
    % Get the indexes for the two leaves
    z_a = Z > Zuq_iter(i);
    z_a_est = Z_est > Zuq_iter(i);
    z_b = ~z_a;
    z_b_est = ~z_a_est;

    % Get optimal bandwidht for both leaves
    bw_a = rdbwselect( Y( z_a ) , X( z_a , 1 ) , 'c', c, 'fuzzy', fuzzy, 'p', optTree.orderPolinomial, ...
        'vce', optTree.varEstimator, 'bwselect', optTree.bwselect );
    bw_b = rdbwselect( Y( z_b ) , X( z_b , 1 ) , 'c', c, 'fuzzy', fuzzy, 'p', optTree.orderPolinomial, ...
        'vce', optTree.varEstimator, 'bwselect', optTree.bwselect );

    % Check if have enough observations
    % above split value
    w_tr_a = feval( kernFnc , ( ( X( z_a & ITT , 1 ) - c ) ./ bw_a.h_r ) ) ./ bw_a.h_r;
    w_co_a = feval( kernFnc , ( ( X( z_a & uITT, 1 ) - c ) ./ bw_a.h_l ) ) ./ bw_a.h_l;
    w_est_tr_a = feval( kernFnc , ( ( X_est(z_a_est & ITT_est, 1 ) - c ) ./ bw_a.h_r ) ) ./ bw_a.h_r;
    w_est_co_a = feval( kernFnc , ( ( X_est(z_a_est & uITT_est, 1 ) - c ) ./ bw_a.h_l ) ) ./ bw_a.h_l;
    % below split value
    w_tr_b = feval( kernFnc , ( ( X( z_b & ITT , 1 ) - c ) ./ bw_b.h_r ) ) ./ bw_b.h_r;
    w_co_b = feval( kernFnc , ( ( X( z_b & uITT, 1 ) - c ) ./ bw_b.h_l ) ) ./ bw_b.h_l;
    w_est_tr_b = feval( kernFnc , ( ( X_est(z_b_est & ITT_est, 1 ) - c ) ./ bw_b.h_r ) ) ./ bw_b.h_r;
    w_est_co_b = feval( kernFnc , ( ( X_est(z_b_est & uITT_est, 1 ) - c ) ./ bw_b.h_l ) ) ./ bw_b.h_l;

    % Get the effective number of observations
    n_tr_a = sum( w_tr_a > 0 );
    n_co_a = sum( w_co_a > 0 );
    n_tr_eff_a = sum( w_est_tr_a > 0 );
    n_co_eff_a = sum( w_est_co_a > 0 );
    n_tr_b = sum( w_tr_b > 0 );
    n_co_b = sum( w_co_b > 0 );
    n_tr_eff_b = sum( w_est_tr_b > 0 );
    n_co_eff_b = sum( w_est_co_b > 0 );

    % If there are enough effective observations
    if ( n_tr_a >= minObs && n_co_a >= minObs && ...
         n_tr_eff_a >= minObs && n_co_eff_a >= minObs && ...
         n_tr_b >= minObs && n_co_b >= minObs && ...
         n_tr_eff_b >= minObs && n_co_eff_b >= minObs )       
    
        % Do the first valid criterion
        %if all( isnan( split_crit( 1:i ) ) )
            sq_a = get_np_split_criterion(Y, X, bw_a, c, z_a, W, cl, optTree, N_tr, N_est, z_a_est );
            sq_b = get_np_split_criterion(Y, X, bw_b, c, z_b, W, cl, optTree, N_tr, N_est, z_b_est );
            split_crit( i ) = sq_a + sq_b;
            % n_tr_last = n_tr_a + n_tr_b;
            % n_co_last = n_co_a + n_co_b;
        % else % Check for bucketing if there are already valid criterion
        %     n_tr_new = ( n_tr_a + n_tr_b ) - n_tr_last;
        %     n_co_new = ( n_co_a + n_co_b ) - n_co_last;
        %     if ( n_tr_new >= optTree.obsBucket && n_co_new >= optTree.obsBucket )
        %         sq_a = get_np_split_criterion(Y, X, bw_a, c, z_a, W, cl, optTree, N_tr, N_est, z_a_est );
        %         sq_b = get_np_split_criterion(Y, X, bw_b, c, z_b, W, cl, optTree, N_tr, N_est, z_b_est );
        %         split_crit( i ) = sq_a + sq_b;
        %         % Update the last valid obs
        %         n_tr_last = n_tr_a + n_tr_b;
        %         n_co_last = n_co_a + n_co_b;
        %     else
        %         exit_flag{ i } = 'Not enough new treatment/control variable.';
        %     end
        % end
    else
        exit_flag{ i } = 'Not enough effective observation!';
    end

end


% Check for possible NaN values
logValid = ~isnan( split_crit );

% If there are any candidates
if any( logValid )
    [ ~ , minID ] = min( split_crit );
    split = Zuq_iter( minID );
else
    split = NaN;
end

end
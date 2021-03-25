%% CT-Athey-Imbens splitting criterion - implementation of C++ code

function split  = split_CT_AI( Y , X , treatment , X_est , tr_est , train_to_est_ratio , optTree )


%% Initialization
n       = numel( Y );
nclass  = 0;
% No weights implemented
wt      = ones( n , 1 );
alpha   = optTree.criteria_weight;
minObs  = optTree.minObs;

%% First calculate the node's effect itself
right_wt = sum( wt );
right_tr = sum( wt .* treatment );
right_sum = sum( Y .* wt );
right_tr_sum = sum( Y .* wt .* treatment );
right_sqr_sum = sum( Y .^ 2 .* wt );
right_tr_sqr_sum = sum( Y.^ 2 .* wt .* treatment );

% Calculate Honest criterion
% temp - treatment effect: difference of means for control and treatment for each split value
temp_1       = right_tr_sum     ./ right_tr - ( right_sum - right_tr_sum ) ./ ( right_wt - right_tr );
% Difference for the square term
tr_var      = right_tr_sqr_sum ./ right_tr - right_tr_sum .^2  ./ ( right_tr .^ 2 );
con_var     = ( right_sqr_sum - right_tr_sqr_sum ) ./ ( right_wt - right_tr ) ...
    - ( right_sum - right_tr_sum ) .* ( right_sum - right_tr_sum ) ./ ...
    ( ( right_wt - right_tr ) .* ( right_wt - right_tr ) );
% Estimate honest criterion for the node
node_effect = alpha .* temp_1.^2 .* right_wt ...
    - ( 1 - alpha ) .* ( 1 + train_to_est_ratio ) .* right_wt ...
    .* ( tr_var ./ right_tr  + con_var ./ ( right_wt - right_tr) );

%% Continuous predictor
if nclass == 0
    % Initialization
    % Number of weights in the left and right parts
    left_wt         = cumsum( wt );
    right_wt        = right_wt - left_wt;
    % Number of treated observation, weighted in the left and right parts
    left_tr         = cumsum( wt .* treatment );
    right_tr        = right_tr - left_tr;
    % Number of observation in the left and right parts
    left_n          = ( 1 : n )';
    right_n         = ( n : -1 : 1 )';
    % Sum of Y
    temp_sY      = Y .* wt;
    left_sum     = cumsum(  temp_sY );
    right_sum    = right_sum - left_sum;
    % Sum of treated means for left and right
    temp_2          = Y .* wt .* treatment;
    left_tr_sum     = cumsum(  temp_2 );
    right_tr_sum    = right_tr_sum - left_tr_sum;
    % square of the mean for the left and right - base
    temp_3        = Y.^2 .* wt;
    left_sqr_sum  = cumsum(  temp_3 );
    right_sqr_sum = right_sqr_sum - left_sqr_sum;
    % square of the mean for the left and right - treatment
    temp_4           = Y.^2 .* wt .* treatment;
    left_tr_sqr_sum  = cumsum(  temp_4 );
    right_tr_sqr_sum = right_tr_sqr_sum - left_tr_sqr_sum;
    
    % Calculate the left part criterion
    left_temp = left_tr_sum ./ left_tr ...
                - ( left_sum - left_tr_sum ) ./ ( left_wt - left_tr );
    left_tr_var =    left_tr_sqr_sum ./ left_tr ...
                   - left_tr_sum.^ 2  ./ left_tr.^2;
    left_con_var =  ( left_sqr_sum - left_tr_sqr_sum ) ./ ( left_wt - left_tr ) ...
                 -  ( left_sum - left_tr_sum ).^2  ...
                    ./ ( left_wt  - left_tr ).^2 ;
    left_effect =  alpha .* left_temp.^2  .* left_wt ...
                - ( 1 - alpha ) .* ( 1 + train_to_est_ratio ) .* left_wt ...
                    .* ( left_tr_var ./ left_tr + left_con_var ./ ( left_wt - left_tr ) );
    
    % Calculate the right part
    right_temp = right_tr_sum ./ right_tr ...
                - ( right_sum - right_tr_sum ) ./ ( right_wt - right_tr );
    right_tr_var =  right_tr_sqr_sum  ./ right_tr ...
                  - right_tr_sum.^ 2  ./ right_tr.^2;
    right_con_var =  ( right_sqr_sum - right_tr_sqr_sum ) ./ ( right_wt - right_tr ) ...
                  -  ( right_sum - right_tr_sum ).^2  ...
                    ./ ( right_wt  - right_tr ).^2 ;
    right_effect =  alpha .* right_temp.^2  .* right_wt ...
                - ( 1 - alpha ) .* ( 1 + train_to_est_ratio ) .* right_wt ...
                .* ( right_tr_var ./ right_tr + right_con_var ./ ( right_wt - right_tr ) );
    
    
    split_values = left_effect + right_effect - node_effect;
    
    %% Select such criterion which satisfies all the criterion
    % Valid splitting criterion
    logNonZero = split_values > 0;
    % Remove nans -> for variance part
    logNaN     = ~isnan( split_values );
    % Summarise valid values
    logValid = logNonZero & logNaN;
    % If there are any candidates
    if any( logValid )
        % Candidate Values
        can_split = split_values( logValid );
        % Now sort them into buckets to smooth irregularities by individual
        % observations -> use buckets
        % Note Athey et al (2016) uses the buckets beforehand. However it does
        % not matter to do it in this way without the for loop.
        % 
        % Find the proper number of observations in each buckets
        can_tr  = treatment( logValid ); 
        sum_tr  = sum( can_tr );
        n_valid = sum( logValid );
        bucketTmp = min( round( sum_tr ./ optTree.obsBucket ) , round( ( n_valid - sum_tr ) ./ optTree.obsBucket ) );
        num_buckets    = max( minObs , min( bucketTmp , optTree.maxBucket ) );
        obs_in_buckets = round( n_valid ./ num_buckets );
        % Runs through the number of observations and select the end of each
        % bucket in the covariate
        log_Bucketing = false( n_valid , 1 );
        n_tr  = 0;
        n_con = 0;
        buc_num_i = NaN( n_valid , 1 );
        ct = 1;
        for i = 1 : n_valid
            % Count the number of treatment or control observation
            if can_tr( i )
                n_tr = n_tr + 1;
            else
                n_con = n_con + 1;
            end
            buc_num_i( i ) = ct;
            % If we have enough observations for both treated and control
            % than we accept as a candidate for the split and restart
            % counting
            if n_tr >= obs_in_buckets && n_con >= obs_in_buckets
                log_Bucketing( i ) = true;
                n_tr = 0;
                n_con = 0;
                ct = ct + 1;
            end
        end
        % Finally check for enough observations in the 'tails'        
        logNumObs = left_tr( logValid ) >= minObs & ...                             Minimum number of treatment in the left side
                    ( left_wt( logValid ) - left_tr( logValid ) ) >= minObs & ...   Minimum number of control in the left side
                    right_tr( logValid ) >= minObs & ...                            Minimum number of treatment in the right side
                    ( right_wt( logValid ) - right_tr( logValid ) ) >= minObs;   %  Minimum number of control in the right side
        % Check in the estimation sample as well
        left_tr_e   = cumsum( tr_est );
        left_con_e  = cumsum( ~tr_est );
        n_est       = numel( tr_est );
        right_tr_e  = n_est - left_tr_e;
        right_con_e = n_est - left_con_e;
        logNum_est  = left_tr_e >= minObs & left_con_e >= minObs & ...
                      right_tr_e >= minObs & right_con_e >= minObs;
        X_e_v       = X_est( logNum_est );
        logEst_valid = X( logValid ) >= X_e_v( 1 ) & X( logValid ) <= X_e_v( end );
        
        
        logFin = log_Bucketing & logNumObs & logEst_valid;
        if any( logFin )
            % Find the best ID
            [ max_crit , max_ID ] = max( can_split( logFin ) );
            % Get the best split value from X
            aux_1 = X( logValid );
            aux_2 = aux_1( logFin );
            % Find out wheter it is treatment or control
            bX = aux_2( max_ID );
            % Find the bucket
            b_bucket = buc_num_i( aux_1 == bX );
            % Splitting value is the largest value in the bucket
            aux_split = aux_1( b_bucket == buc_num_i );
            split = aux_split( end );
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

end
%% Work-horse for updating algorithm:
% Estimate the point estimates and the inverse of Xi'Xi and 

function [ b_y_a , b_y_b , Mi_a , Mi_b , b_w_a , b_w_b ] = wh_sms_est_coeff( ...
                          p , n , n_v , f_obs , l_obs , X0 , Y , W , fuzzy ) 


%% Create variables
Mi_a  = NaN( p + 1 , p + 1 , n_v );
Mi_b  = NaN( p + 1 , p + 1 , n_v );
b_y_a = NaN( n_v , p + 1 );
b_y_b = NaN( n_v , p + 1 );
% If fuzzy it is needed - note the variance for w is not directly needed!
b_w_a = NaN( n_v , p + 1 );
b_w_b = NaN( n_v , p + 1 );




%% Initialize algorithm:

% Below starts from first observation and increase
Mi_b( : , : , 1 ) = inv( 1 ./ f_obs .* X0( 1 : f_obs , : )' * X0( 1 : f_obs , : ) );
MiX_b             =    ( 1 ./ f_obs .* Mi_b( : , : , 1 ) )  * X0( 1 : f_obs , : )';
b_y_b( 1 , : )    = MiX_b * Y( 1 : f_obs , : );

% Above starts from last  observation and decrease
num_obs_a        = n - l_obs + 1;
Mi_a( : , : , end ) = inv( 1 ./ num_obs_a .* X0( l_obs : end , : )' * X0( l_obs : end , : ) );
MiX_a               =    ( 1 ./ num_obs_a .* Mi_a( : , : , end ) )  * X0( l_obs : end , : )';
b_y_a( end , : )    = MiX_a * Y( l_obs : end , : );

% If fuzzy design calculate the treatment equation as well
if fuzzy
    b_w_b( 1 , : )    = MiX_b * W( 1 : f_obs , : );
    b_w_a( end , : )  = MiX_a * W( l_obs : end , : );                            
end
                 

for i = 1 : ( n_v - 1 )
    
    %% Outcome equation
    % update the coefficients and inverse for observation below the split
    [ b_y_b( i + 1 , : ) , Mi_b( : , : , i + 1 ) ] = split_update_OLS( b_y_b( i , : )' ,...
                                                                       Mi_b( : , : , i ) , ...
                                                                       f_obs + i , ...
                                                                       Y( f_obs + i ) , X0( f_obs + i , : ) );

    % update the coefficients and inverse for observation above the split
    [ b_y_a( end - i , : ) , Mi_a( : , : , end - i ) ] = split_update_OLS( b_y_a( end - i + 1 , : )' ,...
                                                                           Mi_a( : , : , end - i + 1 ) , ...
                                                                           num_obs_a + i , ...
                                                                           Y( l_obs - i ) , X0( l_obs - i , : ) );
    if fuzzy
        %% Treatment equation
        %   only beta needs to be updated, inverse is already calculated
        b_w_b( i + 1 , : )   = update_only_b( W( f_obs + i ) , X0( f_obs + i , : ) , ...
                                            b_w_b( i , : )'  , Mi_b( : , : , i + 1 ) , f_obs + i + 1 );
        b_w_a( end - i , : ) = update_only_b( W( l_obs - i ) , X0( l_obs - i , : ) , ...
                                            b_w_a( end - i + 1 , : )' , Mi_a( : , : , end - i ) , num_obs_a + i + 1 );
    end
end



end

function b = update_only_b( Y_new , X_new , b_1 , Mi , n )
    % Update OLS
    b  = b_1 + ( 1 ./ n .* Mi ) * X_new' .* ( Y_new - X_new * b_1 ); 
end
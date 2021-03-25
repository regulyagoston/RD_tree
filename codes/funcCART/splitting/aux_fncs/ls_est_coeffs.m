%% Compute coefficients and M_inverse for split_CT_RDD with simple least squares
%
% Experience: it is faster than lscov if there are more than 10.000 observations or
% less than 3 times the observations


function [ b_y_t_a , b_y_t_b , b_y_c_a , b_y_c_b , ...
           Mi_t_a , Mi_t_b , Mi_c_a , Mi_c_b , ...
           b_w_t_a , b_w_t_b , b_w_c_a , b_w_c_b ] = ls_est_coeffs( ...
                                              p , n_t , n_c , idZ_t , idZ_c , ...
                                              X0_t , X0_c , Y_t , Y_c , W_t , W_c , fuzzy )

%% Compute estimates for valid treatment group
[ b_y_t_a , b_y_t_b , Mi_t_a , Mi_t_b , ...
  b_w_t_a , b_w_t_b ] = est_ls_gen( p , n_t , idZ_t , X0_t , Y_t , W_t , fuzzy );


%% Compute estimates for valid control group
[ b_y_c_a , b_y_c_b , Mi_c_a , Mi_c_b , ...
  b_w_c_a , b_w_c_b ] = est_ls_gen( p , n_c , idZ_c , X0_c , Y_c , W_c , fuzzy );
end

%% Compute estimates for the given group
function [ b_y_a , b_y_b , Mi_a , Mi_b , ...
           b_w_a , b_w_b ] = est_ls_gen( p , n , idZ , X0 , Y , W , fuzzy )

n_v = numel( idZ );
% Create variables:
Mi_a    = NaN( p + 1 , p + 1 , n_v );
Mi_b    = NaN( p + 1 , p + 1 , n_v );
b_y_a   = NaN( n_v , p + 1 );
b_y_b   = NaN( n_v , p + 1 );

% If fuzzy it is needed
b_w_a   = NaN( n_v , p + 1 );
b_w_b   = NaN( n_v , p + 1 );


ct = 1;
for i = idZ'
    
    % Number of observations above the cut-off value for treated units
    n_a = n - i;
    % update the coefficients and inverse for observation below the split
    % Inverse
    Mi_b( : , : , ct ) = inv( 1 ./ i   .* X0( 1 : i , : )'       * X0( 1 : i       , : ) );
    Mi_a( : , : , ct ) = inv( 1 ./ n_a .* X0( i + 1 : end , : )' * X0( i + 1 : end , : ) );
    % Mi * X0
    MiX_b = ( 1 ./ i   .* Mi_b( : , : , ct ) ) * X0( 1 : i , : )';
    MiX_a = ( 1 ./ n_a .* Mi_a( : , : , ct ) ) * X0( i + 1 : end , : )';
    % Coefficients for outcome
    b_y_b( ct , : ) = MiX_b * Y( 1 : i , : );
    b_y_a( ct , : ) = MiX_a * Y( i + 1 : end , : );
    
    if fuzzy
        % Coefficients for getting treated
        b_w_b( ct , : ) = MiX_b * W( 1 : i , : );
        b_w_a( ct , : ) = MiX_a * W( i + 1 : end , : );
    end
    
    % New cycle
    ct = ct + 1;    
end


end
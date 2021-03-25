%% This function update efficiently the coefficient vector and the inverse of X'X
%
% b_1   : former coefficient values
% Mi_1  : former inverse of the matrix X'X 
% n_1   : former number of observations
% Y_new : scalar of new outcome
% X_new : vector of new observations


function [ b , Mi , n ] = split_update_OLS( b_1 , Mi_1 , n_1 , Y_new , X_new )

n = n_1 + 1;
% Sherman-Morrison formula
Mi = n ./ n_1 .* ( Mi_1 - ( Mi_1 * ( X_new' * X_new ) * Mi_1 ) ./ ...
                          ( n_1 + X_new * Mi_1 * X_new' ) );
% Update OLS
b  = b_1 + 1 ./ n * Mi * X_new' .* ( Y_new - X_new * b_1 ); 
    

end
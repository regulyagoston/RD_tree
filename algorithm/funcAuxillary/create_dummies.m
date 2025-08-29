%% Create dummy variables from categorical variable
%
% Input:
%   X - Nx1 vector of categorical variable with K categories
%
% Output:
%   D  - NxK 0-1 variable
%   Kn - 1xK containing the categories
%

function [ D , Kn ] = create_dummies( X )

N  = numel( X );
Kn = unique( X );
nK = numel( Kn );

D = zeros( N , nK );
if isnumeric( X )
    for i = 1 : nK
        D( X == Kn( i ) , i ) = 1;
    end
elseif iscellstr( X )
    for i = 1 : nK
        D( strcmpi( X , Kn{ i } ) , i ) = 1;
    end
elseif ischar( X )
    for i = 1 : nK
        D( strcmpi( X , Kn( i ) ) , i ) = 1;
    end
else
    error('Use numeric, cellstr or char vector for X!')
end

end
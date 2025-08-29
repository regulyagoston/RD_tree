%% RDROBUST VCE

function M = rdrobust_vce( d , s , RX , res , C , whh )

k = size( RX , 2 );
M = zeros( k , k );
n = numel( C );

if isempty( C )
    w = whh;
    if d == 0 
        rRX = res .* RX;
        M = rRX' * rRX;
    else
        for i = 1 : ( 1 + d )
            SS = res( : , i ) .* res;
            for j = 1 : ( 1 + d )
                M = M + ( RX .* ( s( i ) .* s( j ) ) .* SS( : , j ) )' * RX;
            end
        end
    end
else
    clusters = unique( C , 'stable' );
    g = numel( clusters );
    w = whh .* ( ( n - 1 ) ./ ( n - k ) ) .* ( g ./ ( g - 1 ) );
    if d == 0
        for i = 1 : g
           ind = C == clusters( i );
           Xi = RX( ind , : );
           ri = res( ind , : );
           Xiri = Xi' * ri;
           M = M + Xiri * Xiri';
        end
    else
        for i = 1 : g
            ind = C == clusters( i );
            Xi = RX( ind , : );
            ri = res( ind, : );
            for l = 1 : ( d + 1 )
                for j = 1 : ( d + 1 )
                    M = M + ( ( Xi' * ( s( l ) .* ri( : , l ) ) ) * ( Xi' * ( s( j ) .* ri( : , j ) ) )' )';
                end
            end
        end
    end  
end
M = w .* M;
end
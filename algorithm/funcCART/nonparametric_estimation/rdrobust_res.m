%% RDROBUST_RES



function res = rdrobust_res( x , y , T , Z , m , hii , vce , matches , dups , dupsid , d )

n = numel( y );
dT = 0;
dZ = 0;
if ~isempty( T )
    dT = 1;
end
if ~isempty( Z )
    dZ = size( Z , 2 );
end
res = NaN( n , 1 + dT + dZ );

if strcmpi( vce , 'nn' )
    for pos = 1 : n
        rpos = dups( pos ) - dupsid( pos );
        lpos = dupsid( pos ) - 1;
        while ( lpos + rpos ) <  min( matches , n - 1 )
            if ( pos - lpos - 1 ) <= 0
                rpos = rpos + dups( pos + rpos + 1 );
            elseif ( pos + rpos + 1 ) > n
                lpos = lpos + dups( pos - lpos - 1 );
            elseif ( x( pos ) - x( pos - lpos - 1 ) ) > ( x( pos + rpos + 1 ) - x( pos ) )
                rpos = rpos + dups( pos + rpos + 1 );
            elseif ( x( pos ) - x( pos - lpos - 1 ) ) < ( x( pos + rpos + 1 ) - x( pos ) )
                lpos = lpos + dups( pos - lpos - 1 );
            else
                rpos = rpos + dups( pos + rpos + 1 );
                lpos = lpos + dups( pos - lpos - 1 );
            end
        end
        ind_J = max( 0 , ( pos - lpos ) ) : min( n , ( pos + rpos ) );
        y_J = sum( y( ind_J ) ) - y( pos );
        Ji = numel( ind_J ) - 1;
        res( pos , 1 ) = sqrt( Ji ./ ( Ji + 1 ) ) .* ( y( pos ) - y_J ./ Ji );
        if ~isempty( T )
            T_J = sum( T( ind_J ) ) - T( pos );
            res( pos , 2 ) = sqrt( Ji ./ ( Ji + 1 ) ) .* ( T( pos ) - T_J ./ Ji );
        end
        if ~isempty( Z )
            for i = 1 : dZ
                Z_J = sum( Z( ind_J , i ) ) - Z( pos , i );
                res( pos , 1 + dT + i ) = sqrt( Ji ./ ( Ji + 1 ) ) .* ( Z( pos , i ) - Z_J ./ Ji );
            end
        end
    end
else
    if strcmpi( vce , 'simple' ) 
        w = 1;
    elseif strcmpi( vce , 'hce-1' ) 
        w = sqrt( n ./ ( n - d ) );
    elseif strcmpi( vce , 'hce-2' ) 
        w = sqrt( 1./(1-hii));
    elseif strcmpi( vce , 'hce-3' ) 
        w = 1./(1-hii);
    else
        error('No such var-covar metric implemented');
    end
    res( : , 1 ) = w.*( y - m( : , 1 ) );
    if dT == 1
        res( : , 2 ) = w .* ( T - m( : , 2 ) );
    end
    if dZ > 0
        for i = 1 : dZ
            res( : , 1 + dT + i ) = w .* ( Z( : , i ) - m( : , 1 + dT + i ) );
        end
    end
end
                    
                    
end
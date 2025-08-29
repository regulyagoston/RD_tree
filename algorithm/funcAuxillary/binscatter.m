%% Basic bin-scatter



function Yb = binscatter( Y , X , xeval )

neval = numel( xeval );

Yb = NaN( neval , 1 );


for i = 1 : neval - 1
    id_bin = X >= xeval( i ) & X < xeval( i + 1 );
    Yb( i ) = mean( Y( id_bin ) );    
end
Yb( end ) = mean( Y( X >= xeval( end ) ) );

end
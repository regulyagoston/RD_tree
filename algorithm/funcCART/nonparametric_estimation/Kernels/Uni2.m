%% Univariate (rectangular) kernel

function y = Uni2( u )

y = ( 0.5 .* ones( size( u ) ) ) .* ( u > -1 ) .* ( u < 1 );

end






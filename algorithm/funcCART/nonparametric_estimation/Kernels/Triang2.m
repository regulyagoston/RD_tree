%% Univariate triangular kernel

function y = Triang2( u )

y = ( 1 - abs( u ) ).*( abs( u ) <= 1 );

end






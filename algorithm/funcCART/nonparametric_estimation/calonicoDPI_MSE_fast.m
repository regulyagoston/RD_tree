%% Calonico et al. Fast Direct Plug-in MSE for boundary
% Using:
%   Triangular Kernel
%   

function [ hOpt , BiasVar ] = calonicoDPI_MSE_fast( Y , X , c )

N  = numel( X ); 
%range = max( X ) - min( X );
quantil_range = prctile( X , 75 ) - prctile( X , 25 );
h0 = 2.576 * min( std( X ) , quantil_range ./ 1.349 ) * N^(-1/5);

%% Kernel
u = ( X - c )./h0;
K = ( 1 - abs( u ) ) .* ( u > -1 ) .* ( u < 1 );
W = K./h0;
ind = W > 0;
% Reallocate
eY = Y( ind );
eX = X( ind );
eW = W( ind );
eN = sum( ind );

%% Variance
eXc = eX - c;
R = [ ones( eN , 1 ) , eXc ];
RsW = R.*sqrt( eW );
iGamma = inv( RsW' * RsW );
RW = R.*eW;
beta = iGamma * ( RW' * eY );
eps = eY - R * beta;
epsRW = eps.*RW;
Sigma = iGamma * ( epsRW' * epsRW ) * iGamma;
m_var = Sigma( end , end );

%% Bias
Hp = [ 1 , h0 ]';
v1 = RW' * ( eXc ./ h0 ) .^ 4;
%v2 = RW' * ( eXc ./ h0 ) .^ 5;
BC1 = Hp .* ( iGamma * v1 );
BC1 = BC1( end );
%BC2 = Hp .* ( iGamma * v2 );
%BC2 = BC2( end );
%% New estimate for bw based on range
% u2 = ( X - c ) ./ range;
% K2 = ( 1 - abs( u2 ) ) .* ( u2 > -1 ) .* ( u2 < 1 );
% ind2 = K2 > 0;
% eN2 = sum( ind2 );
% eY2 = Y( ind2 );
% eX2 = X( ind2 );
% eW2 = K2( ind2 );
% eX2c = eX2 - c;
R2 = [ R , eXc .^ 2 ];
RsW2 = R2.*sqrt( eW );
Gamma2 = RsW2' * RsW2;
beta2 = Gamma2 \ ( ( R2 .* eW )' * eY );
%% One higher order
%R3 = [ R2 , eX2c.^5 ];
%RsW3 = R3.*sqrt( eW2 );
%Gamma3 = RsW3' * RsW3;
%beta3 = Gamma3 \ ( ( R3.*eW2 )' * eY2 );
%% Final estimate
B1 = BC1 .* beta2( 3 );
%B2 = BC2 .* beta3( 6 );
V = N .* h0 .* m_var;

%% Output
BiasVar = ( V ./ ( 4 * B1.^2 ) ) .^ 0.2;
hOpt = BiasVar .* N.^ (-0.2);


end
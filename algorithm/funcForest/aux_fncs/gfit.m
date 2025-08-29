function [ xvals , g_eta ] = gfit( X , sigma , p , nbin , unif_fraction )

if nargin < 3
    p = 2;
    nbin = 1000;
    unif_fraction = 0.1;
end

xvals = linspace( min( min( X ) - 2 .* std( X ) , 0 ) , max( max( X ) + 2 .* std( X ) , std( X ) ) , nbin )';
binw = xvals( 2 ) - xvals( 1 );

zero_idx = find( xvals <= 0 , 1 , 'last' );
noise_kernel = normpdf( xvals ./ sigma ) .* binw ./ sigma;

if zero_idx > 1
    noise_rotate = noise_kernel( [ zero_idx : nbin , 1 : ( zero_idx - 1 ) ] );
else
    noise_rotate = noise_kernel;
end

XX = NaN( nbin , p );
xvals_nlog = xvals >= 0;
for i = 1 : p
    XX( : , i ) = xvals .^ i .* xvals_nlog;
end

x0 = -1 .* ones( 1 , p );
fun = @( eta )neg_loglik( eta , X , XX , xvals , xvals_nlog , unif_fraction , noise_rotate );
eta_hat = fminsearch( fun , x0 );

g_eta_raw = exp( XX * eta_hat' ) .* xvals_nlog;
g_eta_main = g_eta_raw ./ sum( g_eta_raw );
g_eta = ( 1 - unif_fraction ) .* g_eta_main + ...
        unif_fraction .* xvals_nlog ./ sum( xvals_nlog );


end
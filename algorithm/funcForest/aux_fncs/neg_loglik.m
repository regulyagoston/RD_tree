function nllh = neg_loglik( eta , X , XX , xvals , xvals_nlog , unif_fraction , noise_rotate )

g_eta_raw = exp( XX * eta' ) .* xvals_nlog;

if isinf( sum( g_eta_raw ) ) || ( sum( g_eta_raw ) <= 10^-14 )
    nllh = 1000 .* ( size( X , 1 ) + sum( eta .^ 2 ) );
    return;
end

g_eta_main = g_eta_raw ./ sum( g_eta_raw );
g_eta = ( 1 - unif_fraction ) .* g_eta_main + ...
        unif_fraction .* xvals_nlog ./ sum( xvals_nlog );

% Check this!!
f_eta = ifft( fft( g_eta ) .* fft( noise_rotate ) );
Y_interpolated = interp1q( xvals ,-log( max( f_eta , 0.0000001 ) ) , X );
nllh = sum( Y_interpolated );

end
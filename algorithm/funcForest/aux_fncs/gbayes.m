function out = gbayes( x0 , g_est_x , g_est_g , sigma )

    Kx = normpdf( ( g_est_x - x0 ) ./ sigma );
    post = Kx .* g_est_g;
    post = post ./ sum( post );
    out = sum( post .* g_est_x );

end
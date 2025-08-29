function calib_all = calibrateEB( vars , sigma2 )

if sigma2 <= 0 || min( vars ) == max( vars )
    calib_all = max( max( vars ) , 0 );
    return;
end

sigma = sqrt( sigma2 );
[ eb_prior_x , eb_prior_g ] = gfit( vars , sigma );

if numel( vars ) >= 200 
    % If there are many test points, use interpolation to speed up computations
    calib_x = quantile( vars , 0:0.02:1 )';
    xN = numel( calib_x );
    calib_y = NaN( xN , 1 );
    for j = 1 : xN
        calib_y( j ) = gbayes( calib_x( j ) , eb_prior_x , eb_prior_g , sigma );
    end
    calib_all = interp1q( calib_x , calib_y , vars );
    
else
    vN = numel( vars );
    calib_all = NaN( vN , 1 );
    for j = 1 : vN
        calib_all( j ) = gbayes( vars( j ) , eb_prior_x , eb_prior_g , sigma );
    end
end

end

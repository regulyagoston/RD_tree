function coverage = get_coverage_rf( forest_pred, forest_se, trueCATE , ci_val)



%% Calculate the hit ratio:
% Hit ratio for confidence intervals
lratio = ( 1 - ci_val ) / 2;
CI_xs = norminv( [ lratio , ci_val + lratio ] , 0 , 1 );

hitCI = trueCATE >= ( forest_pred + CI_xs( 1 ) .* forest_se ) & ...
        trueCATE <= ( forest_pred + CI_xs( 2 ) .* forest_se );

% Get the coverage ratio
coverage = mean( hitCI );

end
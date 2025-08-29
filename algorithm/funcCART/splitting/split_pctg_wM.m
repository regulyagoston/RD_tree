%% Split function 1 - enforcing the median with percentiles


function split = split_pctg_wM( Zi , nS )

% If 50 is not part of the linspacing then add it mechanically
if mod( 50 / ( 100 / nS ) , 1 ) ~= 0
    split = prctile( Zi , linspace( 0 , 100 , nS - 1 ) );
    % Add the median to the splitting values
    split = [ prctile( Zi , 50 ) , split ];
else
    split = prctile( Zi , linspace( 0 , 100 , nS ) );
end

end
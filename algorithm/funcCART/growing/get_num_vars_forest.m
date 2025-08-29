%% Set the number of candidate variables at each node when using forest

function K_rf = get_num_vars_forest( num_vars, user_value )

if ~isempty( user_value )
    % In case of user supplied value
    K_rf = min( user_value , num_vars );
else
    % It should be number of features / 3 or if it is less than 5 then
    % use 5. If there are less than 5 features, than use all the number of
    % features
    K_rf = min( max( floor( num_vars ./ 3 ) , 5 ) , num_vars );
end


end
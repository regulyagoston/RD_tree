function [ opt_beta , opt_bw, s_oOS, large_tree, rho ] = cross_validate_bws( obj_sample , obj_optCART )


% Initial point is the optimal bandwidth for the whole population
bw0 = rdbwselect( obj_sample.Y_tr, obj_sample.X_tr , 'c', obj_sample.c, 'fuzzy', obj_sample.W_tr, 'p', obj_optCART.orderPolinomial, ...
        'vce', obj_optCART.varEstimator, 'bwselect', obj_optCART.bwselect );
h_0 = bw0.h_l;
%b_0 = bw0.b_l;
rho = 1;%h_0./b_0;

% Minimum: number of observations within the range is equal to the
% pre-set minimum observations
sort_X = sort(obj_sample.X_est - obj_sample.c);
left_X = sort_X( sort_X < 0 );
if obj_optCART.minObs < numel( left_X )
    left_min = -left_X(end-obj_optCART.minObs+1 );
else
    left_min = left_X( 1 );
end

right_X = sort_X( sort_X > 0 );
if obj_optCART.minObs < numel( right_X )
    right_min = right_X(obj_optCART.minObs );
else
    right_min = right_X( end );
end
h_min = min([(left_min+right_min)*2, 0.5*h_0]);% this is just a rule of thumb to avoide problems in CV
h_max = min( max(-min(left_X)+max(right_X))*1.5, h_0*10);

%b_min = rho .* h_min;
%b_max = rho .* h_max;

opts = optimoptions('fmincon','MaxIter', obj_optCART.num_bw);

[opt_bw,s_oOS] = fmincon( @(x)optim_bw(x, rho, obj_sample, obj_optCART),...
    h_0, [], [], [], [], h_min, h_max, [], opts );

[ ~, opt_beta, large_tree  ] = optim_bw( opt_bw, rho, obj_sample, obj_optCART );



end

function [ crit, beta, obj_tree_i ] = optim_bw( bw, rho, obj_sample, obj_optCART )

obj_optCART_i = copyWithChange(obj_optCART,'bw', bw );
obj_optCART_i = copyWithChange(obj_optCART_i,'bw_b', rho .* bw );
obj_tree_i = growRDDtree( obj_sample, obj_optCART_i );
[ beta , crit ] = cross_validate( obj_tree_i , obj_sample , obj_optCART_i );

end
function [ bw, bw_b ] = bw_grid(obj_sample, obj_optCART)


% Best-guess for overall
bw0 = rdbwselect( obj_sample.Y_tr, obj_sample.X_tr , 'c', obj_sample.c, 'fuzzy', obj_sample.W_tr, 'p', obj_optCART.orderPolinomial, ...
        'vce', obj_optCART.varEstimator, 'bwselect', obj_optCART.bwselect );

bw_b0 = bw0.b_l;
bw0 = bw0.h_l;
rho = 1;%bw0/bw_b0;


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
bw_min = min([(left_min+right_min)*2, 0.5*bw0]);% this is just a rule of thumb to avoide problems in CV

bw_max = min( max(-min(left_X)+max(right_X))*1.5, bw0*2);
num_bw_l = ceil(obj_optCART.num_bw/2);
num_bw_h = obj_optCART.num_bw - num_bw_l + 1;
aux_h = linspace(bw0,bw_max,num_bw_h);
bw = [linspace(bw_min,bw0,num_bw_l),aux_h(2:end)];

bw_b = bw * rho;


end


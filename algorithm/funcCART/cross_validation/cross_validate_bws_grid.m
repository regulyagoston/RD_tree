function [ opt_beta , opt_bw, opt_bw_b, s_oOS, bw, beta, bw_b ] = cross_validate_bws_grid( obj_sample , obj_optCART, bw )


if isempty(bw)
    [ bw, bw_b ] = bw_grid(obj_sample, obj_optCART);
end

num_bw = length( bw );

beta = zeros( 1 , num_bw ).*NaN;
s_oOS = zeros( 1 , num_bw ).*NaN;

if obj_optCART.paralell
    parfor i = 1 : num_bw
        obj_optCART_i = copyWithChange(obj_optCART,'bw', bw(i));
        obj_optCART_i = copyWithChange(obj_optCART_i,'bw_b', bw_b(i));
        try
            obj_tree_i = growRDDtree( obj_sample, obj_optCART_i );
            % Modify bandwidth for cross-validation
            [ beta( i ) , s_oOS( i ) ] = cross_validate( obj_tree_i , obj_sample , obj_optCART_i );
        catch
            beta( i ) = NaN;
            s_oOS( i ) = NaN;
        end
    end
else
    for i = 1 : num_bw
        % Optimal bandwidth for the training sample
        obj_optCART_i = copyWithChange(obj_optCART,'bw', bw(i));
        obj_optCART_i = copyWithChange(obj_optCART_i,'bw_b', bw_b(i));
        try
            obj_tree_i = growRDDtree( obj_sample, obj_optCART_i );
            % Modify bandwidth for cross-validation
            [ beta( i ) , s_oOS( i ) ] = cross_validate( obj_tree_i , obj_sample , obj_optCART_i );
        catch
            beta( i ) = NaN;
            s_oOS( i ) = NaN;
        end
    end
end

% Get the minimum value
[ ~ , mID ] = min( s_oOS );
opt_beta = beta( mID );
opt_bw = bw( mID );
opt_bw_b = bw_b( mID );


end
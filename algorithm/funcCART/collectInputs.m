%% Collect the inputs for criterion functions


function varargout = collectInputs( leaf_nodes , tr , varargin )

% Number of inputs and outputs
nV = numel( varargin );
varargout = cell( nV , 1 );
% Number of leaves in the tree
nL   = numel( leaf_nodes );

for i = 1 : nV
    if isfield( leaf_nodes( 1 ).est_struct , varargin{ i } ) || strcmp( varargin{ i } , 'ordered_tau_i' )
        if contains( varargin{ i } , '_j' )
            if any( strncmp( varargin{ i } , 'Mi' , 2 ) || any( strncmp( varargin{ i } , {'var_y','var_t'} , 5 ) ) ...
                    || strncmp( varargin{ i } , 'sigma2_' , 7 ) )
                [ n1 , n2 ] = size( leaf_nodes( 1 ).est_struct.( varargin{ i } ) );
                varargout{ i } = NaN( n1 , n2 , nL );
                t_mi = true;
            else
                varargout{ i } = NaN( nL , 1 );
                t_mi = false;
            end
            chck = 1;
            type = false;
        elseif contains( varargin{ i } , '_i' )
            type = true;
            if tr
                nObs = numel( leaf_nodes( 1 ).logID_tr );
            else
                nObs = numel( leaf_nodes( 1 ).logID_est );
            end
            varargout{ i } = NaN( nObs , 1 );
            chck = NaN( nObs , nL );
        else
            error('CART:collectInputs:invalidInputType',...
                  'Use such field which contains `_i` or `_j` to indicate what type of output are you after!' )
        end
        if type && strcmp( varargin{ i } , 'tau_i' )
            no_ordering = true;
        elseif type
            no_ordering = false;
            if strcmp( varargin{ i } , 'ordered_tau_i' )
                varargin{ i } = 'tau_i';
            end
        end
        % Collect the data
        ct = 1;
        for l = 1 : nL
            if type && tr
                if no_ordering
                    varargout{ i }( ct : ct + leaf_nodes( l ).n_j_tr - 1 ) = leaf_nodes( l ).est_struct.( varargin{ i } );
                    ct = ct + leaf_nodes( l ).n_j_tr;
                else
                    varargout{ i }( leaf_nodes( l ).logID_tr ) = leaf_nodes( l ).est_struct.( varargin{ i } );
                end
                chck( : , l ) = leaf_nodes( l ).logID_tr;
            elseif type && ~tr
                if no_ordering
                    varargout{ i }( ct : ct + leaf_nodes( l ).n_j_est - 1 ) = leaf_nodes( l ).est_struct.( varargin{ i } );
                    ct = ct + leaf_nodes( l ).n_j_est;
                else
                    varargout{ i }( leaf_nodes( l ).logID_est ) = leaf_nodes( l ).est_struct.( varargin{ i } );
                end
                chck( : , l ) = leaf_nodes( l ).logID_est;
            else
                if t_mi
                    varargout{ i }( : , : , l ) = leaf_nodes( l ).est_struct.( varargin{ i } );
                else
                    varargout{ i }( l , 1 ) = leaf_nodes( l ).est_struct.( varargin{ i } );
                end
            end
        end
        if any( any( any( isnan( varargout{ i } ) ) ) ) || ( ~all( sum( chck , 2 ) == 1 ) )
            %error('There are some problems with the given nodes! There is no complete overlap for the observations!')
        end
    else
        error('CART:collectInputs:invalidInputType',...
            ['There is no such property given for nodeProp as: \n',...
            varargin{ i } ] )
    end
end


end
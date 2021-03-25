%% Categorize features
% Finds the needed statistics for the given feature space

function varargout = categorizeFeatures( obj_tree , Z , varargin )

% Number of inputs and outputs and create slots
nV = numel( varargin );
varargout = cell( nV , 1 );
nZ = size( Z , 1 );
for i = 1 : nV
    varargout{ i } = NaN( nZ , 1 );
end
% Get the leaves
lID = findleaves( obj_tree );
nL  = numel( lID );

% Iterates through each leaf
for i = 1 : nL
    % Search for observations from Z in the given leaf -> get the logicals
    if ~isempty( obj_tree.Node{ lID( i ) }.subsetIdx )
        nC = numel( obj_tree.Node{ lID( i ) }.subsetIdx );
        logId_a = false( [ nZ , nC ] );
        for j = 1 : nC
            if obj_tree.Node{ lID( i ) }.LEq( j )
                logId_a( Z( : , obj_tree.Node{ lID( i ) }.subsetIdx( j ) ) <= obj_tree.Node{ lID( i ) }.splitVal( j ) , j ) = true;
            else
                logId_a( Z( : , obj_tree.Node{ lID( i ) }.subsetIdx( j ) ) >  obj_tree.Node{ lID( i ) }.splitVal( j ) , j ) = true;
            end
        end
        logId = all( logId_a , 2 );
    else
        logId = true( [ nZ , 1 ] );
    end
    % Update all output with the required statistics
    for k = 1 : nV
        try
            if contains( varargin{ k } , '_j' )
                val1 = obj_tree.Node{ lID( i ) }.est_struct.( varargin{ k } );
            elseif contains( varargin{ k } , '_i' )
                aux  = obj_tree.Node{ lID( i ) }.est_struct.( varargin{ k } );
                val1 = aux( 1 );
            else
                error('');
            end
        catch
            error('CART:categorizeFeatures:invalidInputType',...
                  [ varargin{ k } ' is an invalid field! \n',...
                   'Use such field which contains `_i` or `_j` to indicate what type of output are you after!' ] )
        end
        varargout{ k }( logId ) = val1;
    end
end


end
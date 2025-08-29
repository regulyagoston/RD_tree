%% Get the tree with the smallest split criterion
% This is to select feature while growing a tree
%
% Input:
%   tr_tree - tree object
%       contains the original tree
%   parent  - scalar (integer)
%       index for the parent of the candidate nodes
%   c_nodes - Nx2 nodeProp object
%       contains the candidate splits by features
%   data    - sample object
%       contains the data
%   optTree - optCART object
%
% Output:
%   bestNodes - 1x2 nodeProp object
%       contains the best fit split feature
%
% Last modified: 2020/02/20
% Agoston Reguly


function bestNodes = select_feature( tr_tree , parent , c_nodes , data , optTree )

%% First get the tree leaves and remove the parent
lf_IDs = findleaves( tr_tree );
lf_IDs( lf_IDs == parent ) = [];

if ~isempty( lf_IDs ) 
    nodes_T = [ tr_tree.Node{ lf_IDs } ];
else % case of rootnode
    nodes_T = [];
end

% Number of candidate nodes
k = size( c_nodes , 1 );
if k > 1
    % Split criterion for each candidate
    spQ_i = NaN( k , 1 );
    for i = 1 : k
        spQ_i( i ) = optTree.crit_split( c_nodes( i , : ) , data , optTree );
    end
    [ ~ , minId ] = min( spQ_i );
    bestNodes = c_nodes( minId , : );
else
    bestNodes = c_nodes;
end

end
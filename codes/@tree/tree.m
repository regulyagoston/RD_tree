classdef tree
%% TREE  A class implementing a tree data structure.
%
% This class implements a simple tree data structure. Each node can only
% have one parent, and store any kind of data. The root of the tree is a
% privilieged node that has no parents and no siblings.
%
% Nodes are mainly accessed through their index. The index of a node is
% returned when it is added to the tree, and actually corresponds to the
% order of addition.
%
% Basic methods to tarverse and manipulate trees are implemented. Most of
% them take advantage of the ability to create _coordinated_ trees: If
% a tree is duplicated and only the new tree data content is modified
% (i.e., no nodes are added or deleted), then the iteration order and the
% node indices will be the same for the two trees.
%
% Internally, the class simply manage an array referencing the node parent
% indices, and a cell array containing the node data. 

% Jean-Yves Tinevez <tinevez@pasteur.fr> March 2012
    
    properties
        % Hold the data at each node
        Node = { [] };
    end
    properties (SetAccess = private)
        
        % Index of the parent node. The root of the tree as a parent index
        % equal to 0.
        Parent = [ 0 ]; %#ok<NBRAK>
        
    end
    
    methods
        
        % CONSTRUCTOR
        
        function [obj, root_ID] = tree(content, val)
            %% TREE  Construct a new tree
            %
            % t = TREE(another_tree) is the copy-constructor for this
            % class. It returns a new tree where the node order and content
            % is duplicated from the tree argument.
            % 
            % t = TREE(another_tree, 'clear') generate a new copy of the
            % tree, but does not copy the node content. The empty array is
            % put at each node.
            %
            % t = TREE(another_tree, val) generate a new copy of the
            % tree, and set the value of each node of the new tree to be
            % 'val'.
            %
            % t = TREE(root_content) where 'root_content' is not a tree,
            % initialize a new tree with only the root node, and set its
            % content to be 'root_content'.
           
            if nargin < 1
                root_ID = 1;
                return
            end
            
            if isa(content, 'tree')
                % Copy constructor
                obj.Parent = content.Parent;
                if nargin > 1 
                    if strcmpi(val, 'clear')
                        obj.Node = cell(numel(obj.Parent), 1);
                    else
                        cellval = cell(numel(obj.Parent), 1);
                        for i = 1 : numel(obj.Parent)
                            cellval{i} = val;
                        end
                        obj.Node = cellval;
                    end
                else
                    obj.Node = content.Node;
                end
                
            else
                % New object with only root content
                
                obj.Node = { content };
                root_ID = 1;
            end
            
        end
        
        
        % METHODS
        
        function [obj, ID] = addnode(obj, parent, data)
            %% ADDNODE attach a new node to a parent node
            % 
            % tree = tree.ADDNODE(parent_index, data) create a new node
            % with content 'data', and attach it as a child of the node
            % with index 'parent_index'. Return the modified tree.
            % 
            % [ tree ID ] = tree.ADDNODE(...) returns the modified tree and
            % the index of the newly created node.
            
            if parent < 0 || parent > numel(obj.Parent)
                error('MATLAB:tree:addnode', ...
                    'Cannot add to unknown parent with index %d.\n', parent)
            end
            
            if parent == 0
                % Replace the whole tree by overiding the root.
                obj.Node = { data };
                obj.Parent = 0;
                ID = 1;
                return
            end
            
            % Expand the cell by
            obj.Node{ end + 1, 1 } = data;
            
            obj.Parent = [
                obj.Parent
                parent ];
            
            ID = numel(obj.Node);
        
        end
        
        function flag = isleaf(obj, ID)
           %% ISLEAF  Return true if given ID matches a leaf node.
           % A leaf node is a node that has no children.
           if ID < 1 || ID > numel(obj.Parent)
                error('MATLAB:tree:isleaf', ...
                    'No node with ID %d.', ID)
           end
           
           parent = obj.Parent;
           flag = ~any( parent == ID );
           
        end
        
        function IDs = findleaves(obj)
           %% FINDLEAVES  Return the IDs of all the leaves of the tree.
           parents = obj.Parent;
           IDs = (1 : numel(parents)); % All IDs
           IDs = setdiff(IDs, parents); % Remove those which are marked as parent
           
        end
        
        function content = get(obj, ID)
            %% GET  Return the content of the given node ID.
            content = obj.Node{ID};
        end

        function obj = set(obj, ID, content)
            %% SET  Set the content of given node ID and return the modifed tree.
            obj.Node{ID} = content;
        end

        
        function IDs = getchildren(obj, ID)
        %% GETCHILDREN  Return the list of ID of the children of the given node ID.
        % The list is returned as a line vector.
            parent = obj.Parent;
            IDs = find( parent == ID );
            IDs = IDs';
        end
        
        function ID = getparent(obj, ID)
        %% GETPARENT  Return the ID of the parent of the given node.
            if ID < 1 || ID > numel(obj.Parent)
                error('MATLAB:tree:getparent', ...
                    'No node with ID %d.', ID)
            end
            ID = obj.Parent(ID);
        end
        function ID = getparents(obj, IDs)
        %% GETPARENT  Return the ID of the parent of the given node.
            if any( IDs < 1 | IDs > numel( obj.Parent ) )
                logID = IDs < 1 | IDs > numel( obj.Parent );
                error('MATLAB:tree:getparents', ...
                      'No node with ID %d.', IDs( logID ) )
            end
            nID = numel( IDs );
            ID  = NaN( nID , 1 );
            for i = 1 : nID
                ID( i ) = obj.Parent( IDs( i ) );
            end
            
        end
                
        function IDs = getsiblings(obj, ID)
            %% GETSIBLINGS  Return the list of ID of the sliblings of the 
            % given node ID, including itself.
            % The list is returned as a column vector.
            if ID < 1 || ID > numel(obj.Parent)
                error('MATLAB:tree:getsiblings', ...
                    'No node with ID %d.', ID)
            end
            
            if ID == 1 % Special case: the root
                IDs = 1;
                return
            end
            
            parent = obj.Parent(ID);
            IDs = obj.getchildren(parent);
        end
        
        function n = nnodes(obj)
            %% NNODES  Return the number of nodes in the tree. 
            n = numel(obj.Parent);
        end
        
    end

end


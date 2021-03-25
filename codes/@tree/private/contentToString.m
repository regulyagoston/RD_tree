function str = contentToString(content)
% Subfunction to turn any content into a decent string
    if isempty(content)
        % Nothing -> print the void symbol
        str = '?';
    elseif ~ischar(content)
        if numel(content) == 1
            % Scalar stuff we most probably be able to print
            if islogical(content)
                if content
                    str = 'true';
                else
                    str = 'false';
                end
                
            elseif isstruct(content)
                % Struct
                fnames = fieldnames(content);
                if numel(fnames) == 1 
                    % Struct with 1 field -> print it
                    str = [ fnames{1} '->' content.(fnames{1}) ];
                else
                    % Print nbr of fields
                    str = [ 'struct.' num2str(numel(fnames)) '_fields' ];
                end
            elseif isa( content , 'nodeProp' )
                
                if isempty( content.LEq )
                    str = 'root node';
                else
                    if content.LEq( end )
                        signN = '<=';
                    else
                        signN = '>';
                    end
                    str = [ 'z_' , num2str( content.subsetIdx( end ) ) ,...
                                    signN , num2str( content.splitVal( end ) , '%1.2f' ) ];
                end
                % Add the estimated value
                try 
                    str = [ str , ' ' , 'tau_hat = ' , num2str( round( content.est_struct.tau_i( 1 ) , 2 ) , '%1.2f' ) ];
                catch
                    str = [ str , ' ' , 'Y_hat = ' , num2str( round( content.est_struct.Y_hat_i( 1 ) , 2 ) , '%1.2f') ];
                end
                
            elseif iscell(content)
                
                % Cell with one element -> append cell and print element
                str = ['cell:' contentToString(content{1}) ];
                
            else
                % Scalar number -> print it
                str = num2str(content);
            end
            
        else
            % Matrix -> print its size
            dims = size(content);
            str = '<';
            for d = dims
                str = [ str num2str(d) 'x' ]; %#ok<AGROW>
            end
            str(end) = ' ';
            str = [str class(content(1)) '>' ];
        end
    else
        % A string -> print it
        str = content;
    end
end

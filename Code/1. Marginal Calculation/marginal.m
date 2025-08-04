function marginal_function = marginal(variable_vector, variable_name,...
bipartite_matrix, factorization_functions, general_summation, general_multiplication)
%MARGINALS Calculates the Marginal function related to one variable 
%----------------------------------------------------------------------------
%                          Input arguments
%----------------------------------------------------------------------------
%   variable_vectort = contains the information about the
%                      deffinition set of the variables x as ordered
%                      elements inside a vector 
%                     
%                      ex: x in {0,1} then definition set is [0,1]
%
%   variable_index = declares the variable from which we want the
%                    marginalization of the function
%
%                    ex: for variable x1 we input 1
%
%   bipartite_matrix = contains information about the bipartime graph 
%                      which relates the factorization functions with their
%                      respective arguments. This logic matrix has ones
%                      when the functions located in the rows relate to
%                      their arguments located in the columns
%                      
%                      ex:       x1   x2
%                             f1  1    1
%                             f2  0    1  
%  
%   factorization_functions = contains all the functions that factorize the
%                             original function we want to marginalize. It
%                             is a cell array of function hundles
%
%                             ex: f(x1,x2,x3) = f1(x1)f2(x2,x3) so {@f1 @f2}
%
%
%   general_summation, 
%   general_multiplication = the general operations of summation and
%                            multiplication that have the properties
%                            forming a finite field. They are function
%                            handles. 
%
%                            !! Also they need to accept multiple arguments
%
%----------------------------------------------------------------------------

    bm = bipartite_matrix;
    bm_size = size(bipartite_matrix);
    functions_size = bm_size(1);
    variables_size = bm_size(2);

    % The Adjaceny matrix relates all nodes of the graph (varables &
    % functions) with each other. It has one when a pair of nodes is
    % connected. Also this matrix must be square and symmetric
    adjacency_matrix = [zeros(functions_size)  bm; bm', zeros(variables_size)'];
    nodenames_functions = "f" + (1:functions_size);
    nodenames_variables = "x" + (1:variables_size);
    nodenames = [nodenames_functions, nodenames_variables];
    nodenames = convertStringsToChars(nodenames);

    
    % Create the graph 
    factor_graph = graph(adjacency_matrix, nodenames);
    plot(factor_graph)
    marginal_function = message(variable_name, []);
    
    
    % Define the recursive messenger function that propagates 
    % the messages inside the factor graph
    function  message_vector = message(current_node, parent_node)

        % Initialization. Take a cell containing all the names of the 
        % neighbors and remove the parent node
        message_vector = zeros(size(variable_vector));
        neighbors_node = neighbors(factor_graph, current_node);
        parent_index = strcmp(parent_node, neighbors_node);
        neighbors_node(parent_index) = [];
        
        
        % If there aren't any neighbors excluding the parent
        % we are in a leaf node and we terminate the process
        if isempty(neighbors_node)

            if current_node(1) == 'x'
                message_vector = ones(size(variable_vector));
                
            elseif current_node(1) == 'f'
                % If the function itself is leaf node then it must have only one argument
                function_index = str2double(current_node(2 : end));
                message_vector = factorization_functions{function_index}(variable_vector);
            end 
            
            return    
        end
        
        
        % General Recursion. Calculate all message vectors incoming from neighbors
        neighbors_length = length(neighbors_node);
        neighbors_message = cell(1, neighbors_length);
        for index = 1 : neighbors_length
            neighbors_message{index} = message(neighbors_node{index}, current_node);
        end
        
        if current_node(1) == 'x'
            message_vector = arrayfun(general_multiplication, neighbors_message{:});

        elseif current_node(1) == 'f'
            % Numbering of variable, function and argument of the function
            variable_index = str2double(parent_node(2:end));
            function_index = str2double(current_node(2:end));
            argument_index = sum(bipartite_matrix(function_index, 1:variable_index));
            function_kernel = factorization_functions{function_index};
            
            % Create all possible combination of values in two multidimentional grids 
            message_grid = cell(1, neighbors_length);
            variable_grid = cell(1, neighbors_length + 1);
            [message_grid{:}] = ndgrid(neighbors_message{:});
            [variable_grid{:}] = ndgrid(variable_vector);
            
            % Calculate the multiplications and the function onto these grids
            message_product_matrix = arrayfun(general_multiplication, message_grid{:});
            function_kernel_matrix = arrayfun(function_kernel, variable_grid{:});
            
           
            % Put the dimention of the parent variable as last for proper multiplication
            slice_index = repmat({':'}, ndims(function_kernel_matrix), 1);
            dimention_order = [setdiff(1:ndims(function_kernel_matrix), argument_index), argument_index];
            function_kernel_matrix = permute(function_kernel_matrix, dimention_order);
          
            % Calculate the sum over a specified variable value of the parent node
            for index = 1 : length(message_vector)
                slice_index{end} = index;
                values_matrix = arrayfun(general_multiplication, function_kernel_matrix(slice_index{:}), message_product_matrix);
                message_vector(index) = general_summation(values_matrix(:));                
            end
        end
    end
    
end


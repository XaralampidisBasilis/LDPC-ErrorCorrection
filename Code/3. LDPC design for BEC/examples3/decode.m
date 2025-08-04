function [transmited_codeword, codeword_found] = decode(received_codeword, receiver_symbols, Channel_Matrix, ...
Bipartite_Matrix, channel_messages)
%DECODE Error Corects the received codeword in order to extract the
%transmited codeword
%                               Input arguments
%--------------------------------------------------------------------------------
%{              

    received_codeword = the codeword that the receiver detects with errors
                        Y = (y1 y2 y3 ... yn)
                        ex: Y = '01?101?010'
         
    receiver_symbols = a char array contaning the symbols that the receiver
                       has in proper order as the Channel_Matrix
                       ex: y in ['0' '1' '?'] = '01?'
   
    transmiter symbols = it has a default value of '0' of '1'
                       ex: x in ['0' '1'] = '01'            

    Channel_Matrix = a matrix containing the dependet probabilities of the
                     transmition channel, Channel_Matrix(i,j) = PY|X(yj|xi)
                     ex:   xy   0    1    ?    
                           0    1-å  0    å
                           1    0    1-å  å

    Bipartite_Matrix = contains information about the bipartime graph 
                       which relates the factorization functions with their
                       respective arguments. This logic matrix has ones
                       when the functions located in the rows relate to
                       their arguments located in the columns
                       Synonim to parity check matrix 
                      
                       ex:       x1   x2
                              f1  1    1
                              f2  0    1 

%}
%                               Output arguments
%--------------------------------------------------------------------------------
%{
    transmited_codeword = the information codeword that the trasmier sends 
                          through the channel X = (x1 x2 x3 ... xn)
                          ex: X = '0101011010'

    codeword_found =  a logical flag that declares if the codeword
                      was found or we have an error

%}

    connected_fun_to_var = logical(Bipartite_Matrix);
    [function_number, variable_number] = size(connected_fun_to_var);
    function_indices = 1 : function_number;
    variable_indices = 1 : variable_number;
    
    % Compute the channel messages to variables from the Channel_Matrix
    if isempty(channel_messages)
        symbol = num2cell(receiver_symbols);
        channel_message_set = log(Channel_Matrix(1,:) ./ Channel_Matrix(2,:));
        map = containers.Map(symbol, channel_message_set); 
        channel_messages = cell2mat(values(map, num2cell(received_codeword)));
    end
    
    % Forward the channel messages to the variables for the iteration to
    % start
    Variable_Messages = repmat(channel_messages, function_number, 1)';
    Variable_Messages(~connected_fun_to_var') = 0; 
    Function_Messages = zeros(size(connected_fun_to_var));

    % Initialization Variables
    threshold_zero = 1e-5;
    iteration_number = 50;
    codeword_found = false;
    transmited_codeword = received_codeword;
    
    for iteration = 1 : iteration_number
            
        
        % Compute the message from a function node to a variable node 
        for function_node = function_indices
            
            % logical indexing of the connected variable nodes
            connected_variables = connected_fun_to_var(function_node, :);
            
            % Take indexing only from the connected variable nodes
            for variable_node = variable_indices(connected_variables)
                
                % Exclude the variable node we are sending the function message
                input_variables = connected_variables;
                input_variables(variable_node) = false;
                function_message = 2 * atanh(prod(tanh(Variable_Messages(input_variables, function_node) / 2), 1));
                
                % Edge case when function node connects with one variable node                
                if isempty(input_variables)
                    function_message = inf; 
                end
                Function_Messages(function_node, variable_node) = function_message;
            end
        end
        
        
        % Apply the error correction in order to recall the transmited_codeword
        decide_messages = sum(Function_Messages, 1) + channel_messages;
        decide_messages(abs(decide_messages) < threshold_zero) = 0;
        transmited_codeword(decide_messages > 0) = '0';
        transmited_codeword(decide_messages < 0) = '1';
        
        
        % Compute the message from a variable node to a function node 
        for variable_node = variable_indices
            
            % logical indexing of the connected function nodes
            connected_functions = connected_fun_to_var(:, variable_node);
            
            % Take indexing only from the connected function nodes
            for function_node = function_indices(connected_functions)
                
                % Exclude the function node we are sending the variable message
                input_functions = connected_functions;
                input_functions(function_node) = false;
                variable_message = sum(Function_Messages(input_functions, variable_node), 1) + channel_messages(variable_node);
                Variable_Messages(variable_node, function_node) = variable_message;
            end
        end
      
        
        % Terminate the process if a codeword found
        if transmited_codeword <= '1'
            codeword = transmited_codeword - '0';    
            result = ~any(mod(Bipartite_Matrix * codeword', 2));
            if result
                codeword_found = true;
                break;
            end
        end
    end
    
    
end


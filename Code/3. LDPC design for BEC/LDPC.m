function Bipartite_Matrix_Ensemble  = LDPC(variable_degree_histogram, function_degree_histogram, ensemble_number)
%LDPC Creates an Ensemble of LDPC Parity Check Matrices 
%                               Input arguments
%--------------------------------------------------------------------------------
%{     
    Variable_Degree_Distribution = vector containing the number of variable nodes 
                                   Ni with degree i where i is the indexing
                                   of the vector. It can be interpreted as the
                                   distribution of variable node degrees.

    Function_Degree_Distribution = vector containing the number of function nodes 
                                   Ni with degree i where i is the indexing
                                   of the vector. It can be interpreted as the
                                   distribution of function node degrees.

    iterations = a number that declares the size of the ensemble of LDPC
                 parity check matrices

    The Degrees Distribution Pair must satisfy a certein condition. The
    number of edges exiting variable nodes must be the same to the edges 
    entering the function nodes 
%}
%                               Output arguments
%--------------------------------------------------------------------------------
%{
    Bipartite_Matrix_Ensemble = cell containing a number of LDPC Parity Check
                                Matrices aka Bipartite Matrices generated
                                for the algorithm 
%}
    
    % Take the number of variable and function nodes
    variable_number = sum(variable_degree_histogram);
    function_number = sum(function_degree_histogram);
    variable_degree_histogram = variable_degree_histogram(:);
    function_degree_histogram = function_degree_histogram(:);
    
    % Check if edges exiting variables nodes are the same number 
    % with edges entering function nodes else declare an error 
    variable_edges_number = (1 : length(variable_degree_histogram)) * variable_degree_histogram;
    function_edges_number = (1 : length(function_degree_histogram)) * function_degree_histogram;
    if variable_edges_number ~= function_edges_number 
        error('Incorect degree distribution pair, no compatible number of edges');
    end
    
    % Create matrices that will contain information about the 
    % sockets each variable and function node with have
    socket_number = variable_edges_number;    
    Variable_Sockets = zeros(variable_number, socket_number);
    Function_Sockets = zeros(function_number, socket_number);

    % Assign randomly how many degrees each variable and function node will
    % have depending on the input degrees distribution pair
    variable_permutation = randperm(variable_number);
    variable_degree = repelem((1 : length(variable_degree_histogram)), variable_degree_histogram);
    variable_degree(variable_permutation) = variable_degree;
    
    function_permutation = randperm(function_number);
    function_degree = repelem((1 : length(function_degree_histogram)), function_degree_histogram);
    function_degree(function_permutation) = function_degree;
    
    % Assign to each variable and function node sockets depending to the
    % degrees each one has 
    for variable_node = 1 : variable_number
        degree = variable_degree(variable_node);
        variable_socket_indices = (1 : degree) + sum(variable_degree(1 : variable_node - 1));
        Variable_Sockets(variable_node, variable_socket_indices) = 1;
    end
    for function_node = 1 : function_number
        degree = function_degree(function_node);
        function_socket_indices = (1 : degree) + sum(function_degree(1 : function_node - 1));
        Function_Sockets(function_node, function_socket_indices) = 1;
    end
    
    % Create the connections between variable and function
    % sockets randomly with each iteration and compute the ensemble of 
    % LDPC partiry check matrices 
    socket_indices = 1 : socket_number;
    Bipartite_Matrix_Ensemble = cell(ensemble_number, 1);
    
    for iteration = 1 : ensemble_number
        Socket_Connections = zeros(socket_number);
        socket_permutation = randperm(socket_number);
        linear_indices = sub2ind(size(Socket_Connections), socket_indices, socket_permutation);
        Socket_Connections(linear_indices) = 1;
        Bipartite_Matrix_Ensemble{iteration} = mod(Function_Sockets * Socket_Connections * Variable_Sockets', 2);
    end
    
end


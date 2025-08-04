function [variable_degree_histogram, function_degree_histogram, approximate_delta_gap] = ...
        histogramPairLDPC(variables_number, variable_edges_distribution, function_edges_distribution)
%histogramPairLDPC This function takes the propability distributions on
%edges and creates the discrete histograms on node degrees 
%{
%                               Input arguments
%--------------------------------------------------------------------------------
%{
    variables_number = free parameter of the number of variable nodes 
    
    variable_edges_distribution = the optimal distribution that has the number
                                  of variable edges with a specific degree
    
    function_edges_distribution = the optimal distribution that has the number
                                  of function edges with a specific degree
%}
%                               Output arguments
%--------------------------------------------------------------------------------
%{
    variable_degree_histogram = the optimal histogram that has the number
                                of variable nodes with a specific degree

    function_degree_histogram =  the optimal histogram that has the number
                                 of function nodes with a specific degree
%}
%                                  Algorithm
%--------------------------------------------------------------------------------
%{
    The Cutting-Plane Method for Solving Convex Programs 
    https://www.mathworks.com/help/optim/ug/miqp-portfolio-problem-based.html  
    
    Linear Programming and Mixed-Integer Linear Programming
    https://www.mathworks.com/help/optim/linear-programming-and-mixed-integer-linear-programming.html?s_tid=CRUX_lftnav
%}
%}
    
    
    lambda = variable_edges_distribution;
    rho = function_edges_distribution;
    indexes_vector = @(distribution) (1 : length(distribution))';

    % Compute the real node degree histograms from the edges degree distributions
    L_real = variables_number * (lambda ./ indexes_vector(lambda)) / sum(lambda ./ indexes_vector(lambda));
    edges = sum(indexes_vector(L_real) .* L_real);
    R_real = edges * (rho ./ indexes_vector(rho));
    
    % Compute the optimal discrete approximation of variable degree
    % histogram with round while keeping the sum constant
    L_int = smartRound(L_real);
    edges = indexes_vector(L_int)' * L_int;
    
    % Compute the optimal discrete approximation of function degree histogram
    % with the Cutting-Plane Method for Solving Convex Programs 
    R_int = smartRound(R_real);                             % Initial Values
    int_constains = 1 : length(R_int);                      % Integer Values 
    error_coefficients = [-2 * R_real; 1];                  % Constant Vector 
    Equality_Constrains = [indexes_vector(R_int)', 0];      % Constant Number 
    equality_values = edges;                                % of edges 
    Inequality_Constrains = [];
    inequality_values = [];
    functions_number = sum(R_int);
    lb = zeros(1, length(R_real) + 1);                       % Lower Bounds
    ub = [functions_number * ones(1, length(R_real)), inf];  % Upper Bounds
    
    iterations = 10;
    options = optimoptions('intlinprog','Display','none');
    for iteration = 1 : iterations
        
        % New inequality constaints in each iteration in order to
        % approximate nonlinear constaint with taylor approximations
        Inequality_Constrains = [Inequality_Constrains; 2 * R_int', -1];   
        inequality_values = [inequality_values; R_int' * R_int];         
        
        R_int_slack = intlinprog(error_coefficients, int_constains, Inequality_Constrains, inequality_values, ...
                            Equality_Constrains, equality_values, lb, ub, options);
        R_int = round(R_int_slack(1 : end - 1));
        
        % Measurements
        slack_variable = R_int_slack(end);
        quadratic_error = R_int' * R_int - slack_variable;
        R_error = norm(R_real - R_int);
    end
    
    % Compute the LPDC code parameters for the optimal histograms
    score = @(distribution) (1 ./ (1 : length(distribution))) * distribution;
    lambda_new = (indexes_vector(L_int) .* L_int) / sum((indexes_vector(L_int) .* L_int));
    rho_new = (indexes_vector(R_int) .* R_int) / sum((indexes_vector(R_int) .* R_int));
    error_BP = errorBP(lambda_new, rho_new);
    design_rate = 1 - score(rho_new) / score(lambda_new);
    error_Sha = 1 - design_rate;
    approximate_delta_gap = (error_Sha - error_BP) / (1 - error_BP);

    if (indexes_vector(L_int)' * L_int) == (indexes_vector(R_int)' * R_int)
        variable_degree_histogram = L_int;
        function_degree_histogram = R_int;
    else
        error('Distributions do not have Equal Edges')
    end

    
    
end


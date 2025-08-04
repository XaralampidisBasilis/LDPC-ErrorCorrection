function [Variable_Degree_Distribution, Function_Degree_Distribution, error_BP] = parametersLDPC(...
    edges_number, max_variable_degree, max_function_degree, erasure_probability)
%GENERATEDISTRIBUTIONS Summary of this function goes here
%   Detailed explanation goes here


    score = @(distribution) (1 ./ (1 : length(distribution))) * distribution;
    polynomial = @(distribution, x) x .^ (0 : length(distribution) - 1) * distribution;
    variable_coefficients = @(function_polynomial, i)  erasure_probability * ((1 - function_polynomial) .^ (i - 1));
    function_coefficients = @(variable_polynomial, i)  -((1 - erasure_probability * variable_polynomial) .^ (i - 1));

    
    %                                       Initialization
    %________________________________________________________________________________________________
    
    % Free parameters related to coding
    %-----------------------------------
    max_iterations = 100;
    difference_threshold = 1e-5;  % the difference of the design rates in each iteration
    round_decimal = 8;   % the position of decimals we want to round from there after

    % Parameters related to the theory
    %----------------------------------
    x_resolution = 200;     % number of resolution points to approach the interval
    x_step = erasure_probability / x_resolution;
    x_discrete = (0 : x_step : erasure_probability)'; 
       
    % Choose linear programing algorithm 'dual-simplex' 'interior-point' 'interior-point-legacy'
    %--------------------------------------------------------------------------------------------
    options = optimoptions('linprog','Algorithm','dual-simplex');
    
    % Initialize parameters and randomly choose second distribution
    %---------------------------------------------------------------
    design_rate = 0;
    iteration = 0;
    function_distribution = rand(max_function_degree, 1);  
    function_distribution  = function_distribution / sum(function_distribution); 
    
    while true
        disp('running')

        %                            Optimize First Distribution
        %____________________________________________________________________________________
        
        % Optimization Score Coefficients 
        %--------------------------------
        score_vector = -(1 ./ (1 : max_variable_degree));

        % Inequality Conditions 
        %-----------------------
        % Stability condition to avoid problems near x to zero given from
        % the book Modern Coding Theory
        condition_1 = zeros(1, max_variable_degree); 
        condition_1(2) = 1;
        stability_bound = (0 : max_function_degree - 1) * function_distribution;
        stability_bound = 1 / (erasure_probability * stability_bound);
        stability_bound(stability_bound > 1) = 1;

        % this coresponds to the theory f(x,l,r) <= 0 for x in [0 1]
        function_polynomial = polynomial(function_distribution, 1 - x_discrete);
        Function_Conditions_Matrix = variable_coefficients(function_polynomial, 1 : max_variable_degree);
        
        % Combine the conditions
        Inequality_Conditions_Matrix = [Function_Conditions_Matrix; condition_1];
        inequality_values = [x_discrete; stability_bound];

        
        % Equality Conditions
        %---------------------
        condition_2 = ones(1, max_variable_degree);  % The sum must equal to one

        condition_3 = zeros(1, max_variable_degree); % The one degree probability is zero 
        condition_3(1) = 1;

        Equality_Conditions_Matrix = [condition_2; condition_3];
        equality_values = [1; 0];

        % Some pre processing
        %---------------------
        % get rid of some really small numbers that are leftover trash  
        Inequality_Conditions_Matrix = round((10 ^ round_decimal).* Inequality_Conditions_Matrix) ./ (10 ^ round_decimal);
        Inequality_Conditions_Matrix(1, 1) = 0;   % correct the matlab calculation 0^0 = 1 to 0^0 = 0;
     
        
        % Lower and Upper Bounds 
        %------------------------
        lower_bounds = zeros(max_variable_degree, 1);
        upper_bounds = ones(max_variable_degree, 1);

        % Linear Optimization Problem 
        %-----------------------------
        variable_distribution = linprog(score_vector, Inequality_Conditions_Matrix, inequality_values, ...
                                        Equality_Conditions_Matrix, equality_values, lower_bounds, upper_bounds, options);


        %                            Optimize Second Distribution
        %____________________________________________________________________________________________
        
        % Optimization Score Coefficients
        %---------------------------------
        score_vector = (1 ./ (1 : max_function_degree));

        % Inequality Conditions
        %-----------------------
        % this corresponds to the theory g(y,l,r) <= 0 for y in [0 1]
        variable_polynomial = polynomial(variable_distribution, x_discrete);
        Inequality_Conditions_Matrix = function_coefficients(variable_polynomial, 1 : max_function_degree);
        inequality_values = x_discrete - 1;

        % Equality Conditions
        %---------------------
        Equality_Conditions_Matrix = ones(1, max_function_degree); % summation to one
        equality_values = 1;
        
        % Some pre processing
        %---------------------
        % Againg avoid leftover trash 
        Inequality_Conditions_Matrix = round((10 ^ round_decimal).* Inequality_Conditions_Matrix) ./ (10 ^ round_decimal);      
       
        % Lower and Upper Bounds
        %------------------------
        lower_bounds = zeros(max_function_degree, 1);
        upper_bounds = ones(max_function_degree, 1);

        % Linear Optimization Problem
        %-----------------------------
        function_distribution = linprog(score_vector, Inequality_Conditions_Matrix, inequality_values, ...
                                        Equality_Conditions_Matrix, equality_values, lower_bounds, upper_bounds, options);


        %                                 Termination Conditions
        %_____________________________________________________________________________________________
      
        iteration = iteration + 1;
        design_rate_old = design_rate;
        design_rate = 1 - score(function_distribution) / score(variable_distribution);

        if abs(design_rate - design_rate_old) < difference_threshold
            disp('distibitions stabilized'), break
        elseif iteration > max_iterations
            disp('max iterations reached'), break
        end
    end

    %                                 Find the Error Threshold
    %__________________________________________________________________________________________________
    error_BP = errorBP(variable_distribution, function_distribution);

    
    Variable_Degree_Distribution = round(edges_number * (variable_distribution ./ (1 : max_variable_degree)'));
    Function_Degree_Distribution = round(edges_number * (function_distribution ./ (1 : max_function_degree)'));

end


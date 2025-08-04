function [variable_edges_distribution, function_edges_distribution, error_BP, delta_gap] = ...
distributionPairLDPC(max_variable_degree, max_function_degree, erasure_probability)

%distributionPairLDPC Summary of this function goes here
%   Detailed explanation goes here


    score = @(distribution) (1 ./ (1 : length(distribution))) * distribution;
    probability = @(ravg) floor(ravg) * (floor(ravg) + 1 - ravg) / ravg;
    polynomial = @(distribution, x) [ones(size(x)), (x .^ (1 : length(distribution) - 1))] * distribution; % correct the matlab calculation 0^0 = 1 to 0^0 = 0
    variable_coefficients = @(function_polynomial, i) ((1 - function_polynomial) .^ (i - 1));
    
    % Save the input data to the text file results
    results = fopen( 'results.txt', 'wt+' );
    formatSpec = ['Inputs\n\n','max_variable_degree = %u  ', 'max_function_degree = %u  ','erasure_probability = %f\n\n'];
    fprintf(results, formatSpec,max_variable_degree,max_function_degree,erasure_probability);
    
    %                                  Code Parameters
    %________________________________________________________________________________________________
    
    ravg_step = 0.5;
    x_resolution = 100;           
    x_step = erasure_probability / x_resolution;
    x_discrete = (0 : x_step : erasure_probability)'; 
   
    % linprog algorithms 'dual-simplex' 'interior-point' 'interior-point-legacy'
    options = optimoptions('linprog','Algorithm','dual-simplex','Display','none');
    
    % Save the Parameters to text file results
    formatSpec = ['Code Parameters\n\n', 'ravg_step = %f  ', 'x_resolution = %u  ', '\n\nResults\n\n'];
    fprintf(results, formatSpec, ravg_step, x_resolution);
    
    %                                Initialization 
    %_____________________________________________________________________________________________
    max_design_rate = 0;
    
    for ravg = 1 : ravg_step : max_function_degree - 1
        
        % Initialize second distribution with a special form 
        function_edges_distribution = zeros(max_function_degree, 1);
        function_edges_distribution(floor(ravg)) = probability(ravg);
        function_edges_distribution(floor(ravg) + 1) = (1 - probability(ravg));
        
        try
            
            %                            Optimize First Distribution
            %____________________________________________________________________________________

            % Optimization Score Coefficients 
            score_vector = -(1 ./ (1 : max_variable_degree));

            % Inequality Conditions 
            % Stability condition to avoid problems near x to zero given from the book Modern Coding Theory
            condition_1 = zeros(1, max_variable_degree); 
            condition_1(2) = 1;
            stability_bound = 1 / (erasure_probability * (0 : max_function_degree - 1) * function_edges_distribution);

            % this coresponds to the theory f(x,l,r) <= 0 for x in [0 1]
            function_polynomial = polynomial(function_edges_distribution, 1 - x_discrete);
            Function_Conditions_Matrix = variable_coefficients(function_polynomial, 1 : max_variable_degree);

            % Combine the conditions
            Inequality_Conditions_Matrix = [Function_Conditions_Matrix; condition_1];
            inequality_values = [x_discrete ./ erasure_probability; stability_bound];

            % Equality Conditions
            condition_2 = ones(1, max_variable_degree);  % The sum must equal to one
            condition_3 = zeros(1, max_variable_degree); % The one degree probability is zero 
            condition_3(1) = 1;
            Equality_Conditions_Matrix = [condition_2; condition_3];
            equality_values = [1; 0];

            % Lower and Upper Bounds 
            lower_bounds = zeros(max_variable_degree, 1);
            upper_bounds = ones(max_variable_degree, 1);

            % Linear Optimization Problem 
            variable_edges_distribution = linprog(score_vector, Inequality_Conditions_Matrix, inequality_values, ...
                                            Equality_Conditions_Matrix, equality_values, lower_bounds, upper_bounds, options);
                
            % Find optimal solution with respect to maximizing design rate
            design_rate = 1 - score(function_edges_distribution) / score(variable_edges_distribution);            
            if design_rate > max_design_rate
                max_design_rate = design_rate;
                optimal_distribution{1} = variable_edges_distribution;
                optimal_distribution{2} = function_edges_distribution;
            end
            
            % Find the Threshold Error and Multiplicative Gap
            error_BP = errorBP(variable_edges_distribution,  function_edges_distribution);
            error_Sha = 1 - design_rate;
            delta_gap = (error_Sha - error_BP) / (1 - error_BP);
            
            % Save the results to the text file 
            formatSpec = ['ravg = %.1f  ','design_rate = %f  ','error_BP = %f  ','delta_gap == %f\n'];
            fprintf(results,formatSpec,ravg,design_rate,error_BP,delta_gap);
        
        catch
            % Some initial conditions fail and we suppress errors
            condition = ' failed to converge';
            formatSpec = ['ravg = %.1f  ','%s\n'];
            fprintf(results,formatSpec,ravg,condition);
            continue
        end
        
    end
    fclose(results);
    
    % Optimal Distributions
    variable_edges_distribution = optimal_distribution{1};
    function_edges_distribution = optimal_distribution{2};
    
    % Find the Threshold Error and Multiplicative Gap of Optimal Distributions
    error_BP = errorBP(variable_edges_distribution,  function_edges_distribution);
    max_design_rate = 1 - score(function_edges_distribution) / score(variable_edges_distribution);
    error_Sha = 1 - max_design_rate;
    delta_gap = (error_Sha - error_BP) / (1 - error_BP);
end


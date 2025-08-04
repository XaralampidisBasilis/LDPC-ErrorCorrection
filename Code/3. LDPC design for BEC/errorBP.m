function error_BP = errorBP(variable_distribution, function_distribution)
%ERRORBP Summary of this function goes here
%   Detailed explanation goes here

    polynomial = @(distribution, x) [ones(length(x),1) x .^ (1 : length(distribution) - 1)] * distribution(:);
    condition_function = @(x, error) error .* polynomial(variable_distribution, 1 - polynomial(function_distribution, 1 - x)) - x;

    error_interval(1, :) = [0  true];
    error_interval(2, :) = [1  false];
    x_step = 0.01;
    zero_threshold = 1e-7;
        
    for iteration = 1 : 50

        % Find the mean error of the last two errors
        error_BP = mean(error_interval(:, 1));
        x_discrete =  (0 : x_step : error_BP)';  % near zero we have problematic performance
        condition = all(condition_function(x_discrete, error_BP) <= zero_threshold);
        
        % Find the new error condition and take the new interval with 
        % the new error and the old error that has different condition
        index = (error_interval(:, 2) == condition);
        error_interval(index, :) = [error_BP, condition];
    end
    
end


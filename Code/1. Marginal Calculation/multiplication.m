function result = multiplication(varargin)
%MULTIPLICATION Summary of this function goes here
%   Detailed explanation goes here
    
    input_arguments = cell2mat(varargin);
    result = 1;
    for input_index = 1 : length(input_arguments)
        result = result .* input_arguments(input_index);
    end
end


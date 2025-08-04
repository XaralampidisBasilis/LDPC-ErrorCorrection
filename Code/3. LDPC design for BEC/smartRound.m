function x_round = smartRound(x)
%SMARTROUND Summary of this function goes here
%   Detailed explanation goes here

    x_floor = floor(x);
    extra_sum = round(sum(x(:))) - sum(x_floor(:));
    fractional_part = x - x_floor;
    [~, sort_index] = sort(fractional_part);
    sort_index = sort_index(end - extra_sum + 1 : end);
    x_round = x_floor;
    x_round(sort_index) = x_floor(sort_index) + 1;
    
end


function I = simpson_inte(array,pace)
steps = length(array);
edge_sum = sum(array(3:2:steps-2));
mid_sum = sum(array(2:2:steps-1));
I = (1/6)*(2*pace)*(array(1) + 2*edge_sum + 4*mid_sum + array(steps));
end
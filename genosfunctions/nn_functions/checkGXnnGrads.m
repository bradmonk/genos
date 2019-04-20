function [gradnet, gradnum] = checkGXnnGrads(lambda)
%CHECKNNGRADIENTS Creates a small neural network to check the
%backpropagation gradients
%   CHECKNNGRADIENTS(lambda) Creates a small neural network to check the
%   backpropagation gradients, it will output the analytical gradients
%   produced by your backprop code and the numerical gradients (computed
%   using computeNumericalGradient). These two gradient computations should
%   result in very similar values.
%

if ~exist('lambda', 'var') || isempty(lambda)
    lambda = 0;
end

input_layer_size = 3;
hidden_layer_size = 5;
num_labels = 3;
m = 5;

% We generate some 'random' test data
Theta1 = debugGXInitW(hidden_layer_size, input_layer_size);
Theta2 = debugGXInitW(num_labels, hidden_layer_size);

% Reusing debugInitializeWeights to generate X
X  = debugGXInitW(m, input_layer_size - 1);
y  = 1 + mod(1:m, num_labels)';

% Unroll parameters
nn_params = [Theta1(:) ; Theta2(:)];

% Short hand for cost function
costFunc = @(p) GXNNCostF(p, input_layer_size, hidden_layer_size, ...
                               num_labels, X, y, lambda);
[cost, gradnet] = costFunc(nn_params);

gradnum = compGXNumericalGrad(costFunc, nn_params);


disp([gradnet gradnum]);
fprintf(['The above two columns you get should be very similar.\n' ...
         '(Left-NN Gradient, Right-Numerical Gradient)\n\n']);

% Evaluate the norm of the difference between two solutions.  
% If you have a correct implementation, and assuming you used EPSILON = 0.0001 
% in computeNumericalGradient.m, then diff below should be less than 1e-9
diff = norm(gradnum-gradnet)/norm(gradnum+gradnet);

fprintf(['If your backpropagation implementation is correct, then \n' ...
         'the relative difference will be small (less than 1e-9). \n' ...
         '\nRelative Difference: %g\n'], diff);

end

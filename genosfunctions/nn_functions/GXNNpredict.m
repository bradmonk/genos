function [p,a,h] = GXNNpredict(Theta1, Theta2, X)
%PREDICT Predict the label of an input given a trained neural network
%   p = PREDICT(Theta1, Theta2, X) outputs the predicted label of X given the
%   trained weights of a neural network (Theta1, Theta2)

sig    = @(z) 1 ./ (1 + exp(-z) );      % sigmoid activation function
sigrad = @(g) sig(g) .* (1-sig(g));     % sigmoid gradient function
sstf   = @(n) 2 ./ (1 + exp(-2*n)) - 1; % sigmoid symetric transfer func

m = size(X, 1);

h1 = sig([ones(m, 1) X] * Theta1');
h2 = sig([ones(m, 1) h1] * Theta2');

h = h2;

[a, p] = max(h, [], 2);

end

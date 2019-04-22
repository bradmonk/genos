function [p,a,h] = HNpredict(t1, t2, X)

sig    = @(z) 1 ./ (1 + exp(-z) );      % sigmoid activation function
sigrad = @(g) sig(g) .* (1-sig(g));     % sigmoid gradient function

m = size(X, 1);

h1 = sig([ones(m, 1) X] * t1');
h2 = sig([ones(m, 1) h1] * t2');

h = h2;

[a, p] = max(h, [], 2);

end

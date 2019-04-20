function [J, grad] = GenomixCostFun(nn_params, ...
                                   input_layer_size, ...
                                   hidden_layer_size, ...
                                   num_labels, ...
                                   X, y, lambda)
%NNCOSTFUNCTION Implements the neural network cost function for a two layer
%neural network which performs classification
%   [J grad] = NNCOSTFUNCTON(nn_params, hidden_layer_size, num_labels, ...
%   X, y, lambda) computes the cost and gradient of the neural network. The
%   parameters for the neural network are "unrolled" into the vector
%   nn_params and need to be converted back into the weight matrices. 
% 
%   The returned parameter grad should be a "unrolled" vector of the
%   partial derivatives of the neural network.
%

% Reshape nn_params back into the parameters Theta1 and Theta2, the weight matrices
% for our 2 layer neural network
Theta1 = reshape(nn_params(1:hidden_layer_size * (input_layer_size + 1)), ...
                 hidden_layer_size, (input_layer_size + 1));

Theta2 = reshape(nn_params((1 + (hidden_layer_size * (input_layer_size + 1))):end), ...
                 num_labels, (hidden_layer_size + 1));

YY = y;
XX = X;

m = size(X, 1);

sig    = @(z) 1 ./ (1 + exp(-z) );     % sigmoid activation function
sigrad = @(g) sig(g) .* (1-sig(g));     % sigmoid gradient function



% ====================== YOUR CODE HERE ======================
%% Part 1: FEEDFORWARD NEURAL NETWORK
%
%         Feedforward the neural network and return the regularized
%         cost in the variable J.


% PREP y FOR COST FUNCTION EVAL

y  = YY;                  % YY has 7000 instances for each 1:2 labels

yd = eye(num_labels);     % create 10x10 identity matrix for labels

y  = yd(y,:);             % Make y a 5000x10 matrix; cols contain 1 or 0



% MAP FROM LAYER-1 TO LAYER-2

X = XX;

a1 = [ones(m, 1) X];    % a1 is input layer; add bias ones col to X matrix (a1: 5000x401)

z2 = a1 * Theta1';      % Coverts to matrix of 5000 examples x 25 thetas (z2: 5000x25)

a2 = sig(z2);       % Sigmoid function converts to p between 0 to 1 (a2: 5000x25)



% MAP FROM LAYER 2 TO LAYER 3

a2 = [ones(m, 1) a2];     % Add ones to the h1 data matrix (a2: 5000x26)

z3 = a2*Theta2';          % Converts to matrix of 5000 exampls x num_labels (z3: 5000x2)

a3 = sig(z3);         % Sigmoid function converts to p between 0 to 1 (a3: 5000x2)




% EVALUATE UNREGULARIZED LOGISTIC COST FUNCTION

logisf = y .* log(a3) + (1-y) .* log(1-a3);   % y as a matrix uses element-wise product (.*)

Jun = (-1/m) .* sum(sum(logisf));             % Cost function without regularization



% REGULARIZE COST FUNCTION

T1 = Theta1(:,2:end);    % 25x400

T2 = Theta2(:,2:end);    % 10x25

J = Jun + (lambda/(2*m)) .* (sum(sum(T1.^2)) + sum(sum(T2.^2)));






%% Part 2: BACKPROPAGATION
%
%         Implement the backpropagation algorithm to compute the gradients
%         Theta1_grad and Theta2_grad. You should return the partial derivatives of
%         the cost function with respect to Theta1 and Theta2 in Theta1_grad and
%         Theta2_grad, respectively.
%
%         Note: The vector y passed into the function is a vector of labels
%           containing values from 1..K. You need to map this vector into a 
%           binary vector of 1's and 0's to be used with the neural network
%           cost function. We recommend implementing backpropagation using a 
%           for-loop over the training examples if you are implementing it 
%           for the first time.



D3 = a3 - y;

D2 = D3 * T2 .* sigrad( z2 );

Delta2  = (1/m) .* (D3' * a2); % Same size as Theta2_grad (10x26)

Delta1  = (1/m) .* (D2' * a1); % Same size as Theta1_grad (25x401)


%% Part 3 
%         Implement regularization with the cost function and gradients.
%         Note: You can implement this around the code for backpropagation.
%           That is, compute the gradients for the regularization separately
%           then add them to Theta1_grad and Theta2_grad from Part 2.
%



% =========================================================================

grad = [Delta1(:) ; Delta2(:)];    % Unroll gradients

end

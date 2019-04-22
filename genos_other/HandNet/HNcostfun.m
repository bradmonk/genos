function [J, grad] = HNcostfun(Tin, L1neurons, L2neurons, nLabs, X, y, lambda)


% Reshape np to param Mx Theta1 & Theta2 weights for 2-layer NN

Theta1 = reshape(Tin(1:L2neurons * (L1neurons + 1)), L2neurons, (L1neurons + 1));
Theta2 = reshape(Tin((1 + (L2neurons * (L1neurons + 1))):end), nLabs, (L2neurons + 1));


YY = y;
XX = X;

m = size(X, 1);


%% Part 1: FEEDFORWARD NEURAL NETWORK

% PREP y FOR COST FUNCTION EVAL

y  = YY;             % YY are label instances for each 1:2 labels

yd = eye(nLabs);     % create identity matrix for labels

y  = yd(y,:);        % Make y m-by-n matrix; cols contain 1 or 0



% MAP FROM LAYER-1 TO LAYER-2

X = XX;

a1 = [ones(m, 1) X];   % a1 is input layer; add bias ones col

z2 = a1 * Theta1';     % covert to matrix of m examples x n thetas

a2 = HNsig(z2);        % sigmoid function



% MAP FROM LAYER 2 TO LAYER 3

a2 = [ones(m, 1) a2];  % Add bias

z3 = a2 * Theta2';     % covert to matrix of examples x num_labels

a3 = HNsig(z3);        % sigmoid function n x 2




% EVALUATE UNREGULARIZED LOGISTIC COST FUNCTION

logisf = y .* log(a3) + (1-y) .* log(1-a3);   % y as a matrix uses element-wise product (.*)

Jun = (-1/m) .* sum(sum(logisf));             % cost function without regularization



% REGULARIZE COST FUNCTION

T1 = Theta1(:,2:end);

T2 = Theta2(:,2:end);

J = Jun + (lambda/(2*m)) .* (sum(sum(T1.^2)) + sum(sum(T2.^2)));






%% Part 2: BACKPROPAGATION

D3 = a3 - y;

D2 = D3 * T2 .* HNsigrad( z2 );

Delta2  = (1/m) .* (D3' * a2); % Same size as Theta2_grad

Delta1  = (1/m) .* (D2' * a1); % Same size as Theta1_grad



grad = [Delta1(:) ; Delta2(:)];    % Unroll gradients

end

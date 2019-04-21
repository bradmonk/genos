function [J, grad] = NNCostF(Tin, L1neurons, L2neurons, nLabs, X, y, lambda)


% Reshape np back into param Mx Theta1 & Theta2 weights for 2-layer NN

Theta1 = reshape(Tin(1:L2neurons * (L1neurons + 1)), L2neurons, (L1neurons + 1));
Theta2 = reshape(Tin((1 + (L2neurons * (L1neurons + 1))):end), nLabs, (L2neurons + 1));
% Theta1 = reshape(Tin(1:L2neurons*L1neurons), L2neurons, L1neurons);
% Theta2 = reshape( Tin(L2neurons*L1neurons+1:end), nLabs, L2neurons );


YY = y;
XX = X;

m = size(X, 1);


%% Part 1: FEEDFORWARD NEURAL NETWORK

% PREP y FOR COST FUNCTION EVAL

y  = YY;             % YY has 7000 instances for each 1:2 labels

yd = eye(nLabs);     % create 10x10 identity matrix for labels

y  = yd(y,:);        % Make y a 5000x10 matrix; cols contain 1 or 0



% MAP FROM LAYER-1 TO LAYER-2

X = XX;

a1 = [ones(m, 1) X];   % a1 is input layer; add bias ones col to X matrix (a1: 5000x401)

z2 = a1 * Theta1';     % Coverts to matrix of 5000 examples x 25 thetas (z2: 5000x25)

a2 = SIG(z2);          % Sigmoid function converts to p between 0 to 1 (a2: 5000x25)



% MAP FROM LAYER 2 TO LAYER 3

a2 = [ones(m, 1) a2];  % Add ones to the h1 data matrix (a2: 5000x26)

z3 = a2 * Theta2';     % Converts to matrix of 5000 exampls x num_labels (z3: 5000x2)

a3 = SIG(z3);          % Sigmoid function converts to p between 0 to 1 (a3: 5000x2)




% EVALUATE UNREGULARIZED LOGISTIC COST FUNCTION

logisf = y .* log(a3) + (1-y) .* log(1-a3);   % y as a matrix uses element-wise product (.*)

Jun = (-1/m) .* sum(sum(logisf));             % Cost function without regularization



% REGULARIZE COST FUNCTION

T1 = Theta1(:,2:end);    % 25x400

T2 = Theta2(:,2:end);    % 10x25

J = Jun + (lambda/(2*m)) .* (sum(sum(T1.^2)) + sum(sum(T2.^2)));






%% Part 2: BACKPROPAGATION

D3 = a3 - y;

D2 = D3 * T2 .* SIGRAD( z2 );

Delta2  = (1/m) .* (D3' * a2); % Same size as Theta2_grad (10x26)

Delta1  = (1/m) .* (D2' * a1); % Same size as Theta1_grad (25x401)





grad = [Delta1(:) ; Delta2(:)];    % Unroll gradients

end

%%  PERFORM MACHINE LEARNING USING HOMEBREW METHODS
% LOGISTIC REGRESSION NEURAL NET WITH GRADIENT DESCENT
%
%
% Assume *TRNN* matrix is formatted as follows
%
% ROW-1: AMX.VID of variant for TRNN(1,4:end)
% COL-1: participant ID
% COL-2: participant AD CASE\CTRL status
% COL-3: bias input node 'ones'
%
% TRNN(2:end,3:end) ... the bias + the variant data
% TRNN(2:end,4:end) ... just the variant data



% DEAL OUT VARIANTS AND LABELS FOR TRAINING DATA
ADNN    =  TRNN;
TRAINX  =  ADNN(2:end,4:end);  
TRAINL  =  ADNN(2:end,2) + 1;


% DEAL OUT VARIANTS AND LABELS FOR TESTING DATA
TESTX   =  TENN(2:end,4:end);
TESTL   =  TENN(2:end,2) + 1;



% RANDOMIZE ROWS
randp        = randperm(size(TRAINX,1));
TRAINX       = TRAINX(randp,:);
TRAINL       = TRAINL(randp,:);


% SET NUMBER OF L1 NEURONS EQUAL TO NUMBER OF VARIANTS
ILayerN  = size(TRAINX,2);
nLabels = 2;



% NEURAL NET PARAMETERS
%----------------------------------------------------------------------
lambda = 0.005;         % .001 - .01 is a good operational window
epsInit = 0.22;         % random initial theta weights
maxIters  = 50;         % 20-50 iterations should be sufficient
HLayerN = 35;           % 10 - 50 neurons should be sufficient
%----------------------------------------------------------------------



% ESTABLISH NEURON ACTIVATION FUNCTION (choose 1)
sig    = @(z) 1 ./ (1 + exp(-z) );     % sigmoid activation function
sigrad = @(g) sig(g) .* (1-sig(g));    % sigmoid gradient function
sstf = @(n) 2 ./ (1 + exp(-2*n)) - 1;  % sigmoid symetric transfer func



% PREP NEURAL NET SYNAPTIC WEIGHTS (THETAS) & GRADIENT DECENT COST FUNC
%----------------------------------------------------------------------

% INITIALIZE RANDOM THETA WEIGHTS
initTheta1 = rand(HLayerN, ILayerN+1) * 2 * epsInit - epsInit;
initTheta2 = rand(nLabels, HLayerN+1) * 2 * epsInit - epsInit;


% UNROLL THETA MATRICES
initial_Thetas = [initTheta1(:) ; initTheta2(:)];
nn_Thetas = initial_Thetas;


% ESTABLISH COST FUNCTION
J = HNcostfun(nn_Thetas, ILayerN, HLayerN, nLabels, TRAINX, TRAINL, lambda);



% TRAIN NEURAL NETWORK USING TRAINING DATASET
%----------------------------------------------------------------------
options = optimset('MaxIter', maxIters);

GcostFun = @(T) HNcostfun(T, ILayerN, HLayerN, nLabels, TRAINX, TRAINL, lambda);

[nn_Thetas, cost] = HNfmincg(GcostFun, initial_Thetas, options);
%----------------------------------------------------------------------




% REROLL THETA WEIGHTS
Theta1 = reshape(nn_Thetas(1:HLayerN * (ILayerN + 1)), HLayerN, (ILayerN + 1));

Theta2 = reshape(nn_Thetas((1 + (HLayerN * (ILayerN + 1))):end), nLabels, (HLayerN + 1));



% EVALUATE PERFORMANCE ON **TRAINING** SET
%----------------------------------------------------------------------

[p , a , h] = HNpredict(Theta1, Theta2, TRAINX);

% p : prediction label {1,2}
% a : confidence level of p
% h : confidence level of p and min label



TRAINPCTCORRECT = mean(p == TRAINL);

disp('Percent accuracy on training data:')
disp( TRAINPCTCORRECT )




% EVALUATE NEURAL NETWORK PERFORMANCE ON **TEST** SET
%----------------------------------------------------------------------

[p , a , h] = HNpredict(Theta1, Theta2, TESTX);


TESTPCTCORRECT = mean(p == TESTL);
disp('Percent accuracy on testing data:')
disp( TESTPCTCORRECT )


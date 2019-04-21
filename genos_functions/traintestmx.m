function [TRAINMX,TRAINL,TESTMX,TESTL,varargout] = traintestmx(TRNN, TENN)
%% makeADTRAINTESTMX.m



% STEP-1: Scramble up TRNN and TENN

NNMX = TRNN(2:end,1:end);
N    = size(NNMX,1);           % Total number of people
xi   = randperm(N)';           % Get N random integers in range 1:N
NNMX = NNMX(xi,:);             % Scramble NN matrix


NNMXTe = TENN(2:end,1:end);
N      = size(NNMXTe,1);       % Total number of people
xi     = randperm(N)';         % Get N random integers in range 1:N
NNMXTe = NNMXTe(xi,:);         % Scramble NN matrix






% STEP-2: Separate into training and testing matrices

TRAIN = NNMX;         % Make training matrix

TEST  = NNMXTe;       % Make testing matrix







% STEP-3: MAKE CLASS LABEL MATRIX Nx2

trnLabels = TRAIN(1:end,2);       % Get AD labels from TRAIN mx

SRR.train  = TRAIN(1:end,1);

TR = zeros(size(trnLabels,1),2);  % Make 2-D class label matrix

TR(:,1) = trnLabels==1;           % In col-1 set CASEs index 1's
TR(:,2) = trnLabels==0;           % In col-2 set CTRLs index 1's



tesLabels = TEST(1:end,2);        % Get AD labels from TEST mx

SRR.test  = TEST(1:end,1);


TE = zeros(size(tesLabels,1),2);  % Make 2-D class label matrix

TE(:,1) = tesLabels==1;           % In col-1 set CASEs index 1's
TE(:,2) = tesLabels==0;           % In col-2 set CTRLs index 1's




% STEP-4: STRIP CASE CTRL INFO FROM TRAINING/TESTING MATRICES

TRAINL = TR;               % Set training label matrix to output variable

TESTL  = TE;               % Set testing label matrix to output variable

TRAINMX = TRAIN(:,4:end);  % Strip case/ctrl info from TRAIN and prep output

TESTMX  = TEST(:,4:end);   % Strip case/ctrl info from TEST and prep output


fprintf('  %-4.0f  CASE labels in training dataset \n',sum(TRAINL(:,1)))
fprintf('  %-4.0f  CTRL labels in training dataset \n',sum(TRAINL(:,2)))
fprintf('  %-4.0f  CASE labels in testing  dataset \n',sum(TESTL(:,1)))
fprintf('  %-4.0f  CTRL labels in testing  dataset \n',sum(TESTL(:,2)))


varargout = {SRR};


end
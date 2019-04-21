function [PCAMX,PCASRR] = makePCAMX(ADNN)
%% makePCAMX.m

% STEP-1: Scramble up the entire NNmx

NNMX = ADNN(2:end,1:end);

N = size(NNMX,1);               % Total number of people

xi = randperm(N)';              % Get N random integers in range 1:N

NNMX = NNMX(xi,:);              % Scramble NN matrix

PCAX = NNMX;                    % Make PCA matrix


% STEP-2: MAKE LABEL MATRIX AND PCA MATRIX

PCASRR = PCAX(1:end,1:2);       % Get AD labels from PCAX

PCAMX  = PCAX(:,4:end);         % Strip case/ctrl info from PCAX


end
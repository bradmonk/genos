function [NNTr,NNTe] = randnnmx(tR,COHORTS)
%% makeSUBMX.m



% STEP-1: Scramble up the entire NNmx

NNMX = COHORTS;

N = size(NNMX,1);               % Total number of people

xi = randperm(N)';              % Get N random integers in range 1:N

NNMX = NNMX(xi,:);              % Scramble CC matrix



% STEP-2: Separate into training and testing matrices

Nt = round(N * tR);          % Use Nt people in the training set

NNTr = NNMX(1:Nt,:);         % Make training matrix
NNTe = NNMX(Nt+1:end,:);     % Make testing matrix


end
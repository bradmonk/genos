%% IMPORT TRAINING DATA
clc; close all; clear

cd(fileparts(which('ML_CLASSIFIER.m')));

load('NNSAMPLEDATA.mat')

%%




% [a,b,c] = pca(ADNN(:,4:end));
% 
% scatter(a(:,1),a(:,2))





%% PREP DATA FOR NN CLASSIFIER TRAINING

disp(ADNN(1:8 , 1:6))
% ADNN
%     1       2       3      4      5      6 ...
%------------------------------------------------
%1|   0    CaCo    bias  SNPid  SNPid  SNPid
%2| SRR     1/0       1    1/0    1/0    1/0
%3| SRR     1/0       1    1/0    1/0    1/0
%4| SRR     1/0       1    1/0    1/0    1/0


TRAINMX      = ADNN(2:end,3:end);   % col 3 of ADNN is bias
TRAINYN      = ADNN(2:end,2);

TESTMX       = ADNNt(2:end,3:end);  % col 3 of ADNN is bias
TESTYN       = ADNNt(2:end,2);


% RANDOMIZE SAMPLE ROWS
rtrain         = randperm(size(TRAINMX,1));
TRAINSAMPLES   = TRAINMX(rtrain,:);
TRAINLABELS    = TRAINYN(rtrain,:);

rtest          = randperm(size(TESTMX,1));
TESTSAMPLES    = TESTMX(rtest,:);
TESTLABELS     = TESTYN(rtest,:);


clearvars -except TRAINSAMPLES TRAINLABELS TESTSAMPLES TESTLABELS

%% TRAIN NEURAL NETWORK CLASSIFIER

% nprtool


%% EVALUATE PERFORMANCE ON TRAINED NET USING "TEST" DATASET

% GuessConfidence = net(TESTSAMPLES');
% 
% yn = TESTLABELS' == round(GuessConfidence);
% 
% PctCorrect = mean(yn) * 100;
% 
% fprintf('\nPercent Correct: %2.2f \n', PctCorrect)


%% EXAMINE & REPORT ACCURACY FOR LOW/MED/HIGH CONFIDENCE DATA

LowConfThresh = .15;
MidConfThresh = .27;
HiConfThresh  = .35;

a = abs(GuessConfidence - .5);


LowConf      = a >= LowConfThresh & a < MidConfThresh;
LowConfPCT   = mean(yn(LowConf)) * 100;
LowConfNp    = numel(yn(LowConf));
LowConfPCTp  = (numel(yn(LowConf)) / numel(yn))*100;


MidConf      = a >= MidConfThresh & a < HiConfThresh;
MidConfPCT   = mean(yn(MidConf)) * 100;
MidConfNp    = numel(yn(MidConf));
MidConfPCTp  = (numel(yn(MidConf)) / numel(yn))*100;

HiConf      = a >= HiConfThresh;
HiConfPCT   = mean(yn(HiConf)) * 100;
HiConfNp    = numel(yn(HiConf));
HiConfPCTp  = (numel(yn(HiConf)) / numel(yn))*100;


fprintf('\nPercent accuracy on mild conf (a > %.2f) testing data: >>> %.1f%% <<<\n',...
    LowConfThresh, LowConfPCT)
fprintf('    (includes %.0f of %.0f [%.1f%%] HighConf:Total patients) \n\n',...
    LowConfNp , numel(yn) , LowConfPCTp )

fprintf('\nPercent accuracy on mid conf (a > %.2f) testing data: >>> %.1f%% <<<\n',...
    MidConfThresh, MidConfPCT)
fprintf('    (includes %.0f of %.0f [%.1f%%] HighConf:Total patients) \n\n',...
    MidConfNp , numel(yn) , MidConfPCTp )

fprintf('\nPercent accuracy on high conf (a > %.2f) testing data: >>> %.1f%% <<<\n',...
    HiConfThresh, HiConfPCT)
fprintf('    (includes %.0f of %.0f [%.1f%%] HighConf:Total patients) \n\n',...
    HiConfNp , numel(yn) , HiConfPCTp )


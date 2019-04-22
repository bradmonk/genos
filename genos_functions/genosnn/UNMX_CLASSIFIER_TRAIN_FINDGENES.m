%% GENOMICS CLASSIFIER USING LOGISTIC FUNCTION WITH GRADIENT DESCENT
% Train a classifier to identify Alzheimer's disease
% using SNP profiles.
%% IMPORT TRAINING DATA
clc; close all; clear
% system('sudo purge')

cd(fileparts(which('UNMX_CLASSIFIER_TRAIN_FINDGENES.m')));

load('ADNETS2.mat')
%load('UNNETS.mat')


disp('ADNN:')
disp(ADNN(1:8 , 1:6))

% [xlsN] = xlsread('TopG.xlsx');



%% PREP DATA FOR NN CLASSIFIER TRAINING
clearvars -except ADNN ADNNt ADGENES ADMX xlsN


for nn = 1:10


% SETABLE PARAMETERS
% if nn < 11
% %-----------------------------------
% UseNgenes = 70;   % 140 genes work best
% maxiters  = 21;    % 21 iterations is perfect
% lambda = 0.006;    % .004 - .008 is best window
% L2neurons = 23;    % 17 - 30 neurons works best
% useBias = 1;       % always use bias
%-----------------------------------
% else
% RANDOMLY GENERATE PARAMETERS
%-----------------------------------
% UseNgenes = randi(50,1) + 45;
% maxiters  = randi(15,1) + 20;
% lambda = randperm(10,1)./1000;
% L2neurons = randi(20)+15;
% useBias = 1;
%-----------------------------------
% end

xx = [...
67	21	30	0.008;
82	21	17	0.004;
77	21	23	0.006];

UseNgenes = xx(mod(nn,3)+1,1);
maxiters  = xx(mod(nn,3)+1,2);
L2neurons = xx(mod(nn,3)+1,3);
lambda    = xx(mod(nn,3)+1,4);
useBias = 1;



INPUTMX       = ADNN(2:end,3:end);  % col 3 of ADNN is bias
INPUTYN       = ADNN(2:end,2);

ADNN(1:9,1:9)
INPUTMX(1:9,1:9)

% size(INPUTMX)
% size(INPUTYN)

% ADNN
%     1       2       3      4      5      6 ...
%------------------------------------------------
%1|   0       0       0  SNPid  SNPid  SNPid
%2| SRR    CaCo    bias    1/0    1/0    1/0
%3| SRR    CaCo    bias    1/0    1/0    1/0
%4| SRR    CaCo    bias    1/0    1/0    1/0


% INPUTMX
%      1      2      3      4
%------------------------------------------------
%1| bias    1/0    1/0    1/0
%2| bias    1/0    1/0    1/0
%3| bias    1/0    1/0    1/0

% xlsG = xlsN-2;

INPUTMXdims   = size(INPUTMX);
INPUTYNdims   = size(INPUTYN);

% RANDOMIZE ROWS
randp         = randperm(INPUTMXdims(1));
INPUTMX       = INPUTMX(randp,:);
INPUTYN       = INPUTYN(randp,:);


disp('INPUTMX:')
disp(INPUTMX(1:9,1:8));    % now col 1 of INPUTMX is bias

% y = NNtoolOUT2(INPUTMX);
% mean(INPUTYN == round(y))


clearvars -except ...
nn NNET ADNN ADNNt ADGENES INPUTMX INPUTYN INPUTMXbias ...
UseNgenes UseGenes maxiters lambda L2neurons useBias xlsN xlsG ...
initTheta1 initTheta2;



%% LOGISTIC REGRESSION NEURAL NET WITH GRADIENT DECENT

sig    = @(z) 1 ./ (1 + exp(-z) );     % sigmoid activation function
sigrad = @(g) sig(g) .* (1-sig(g));     % sigmoid gradient function


% UseGenes = 1:size(INPUTMX,2);
% UseGenes = [1 xlsG'];
% UseGenes = [1 2 randperm( size(INPUTMX,2)-2, UseNgenes)+2 ];
% UseGenes = [1 xlsG' randperm( size(INPUTMX,2)-1, UseNgenes)+1 ];
UseGenes = [1:65 (randperm(size(INPUTMX,2)-65,UseNgenes-65)+UseNgenes-65)];

X  = INPUTMX( : , UseGenes );         % INCLUDE BIAS COL

y = INPUTYN + 1;    % ADD 1 SO LABELS ARE {1,2}


m = size(X,1);      % ROWS: NUMBER OF TRAINING EXAMPLES (PEOPLE)
n = size(X,2);      % COLS: NUMBER OF FEATURES (SNPs)


L1neurons  = n;
num_labels = 2;

% keyboard

% RANDOMIZE INITIAL THETA WEIGHTS
eps_init = 0.12;
% if nn==1

initTheta1 = rand(L2neurons, L1neurons+1) * 2 * eps_init - eps_init;
initTheta2 = rand(num_labels, L2neurons+1) * 2 * eps_init - eps_init;

% initTheta1 = rand(L2neurons, L1neurons) * 2 * eps_init - eps_init;
% initTheta2 = rand(num_labels, L2neurons) * 2 * eps_init - eps_init;


% end

% UNROLL THETA PARAMETERS 
initial_Thetas = [initTheta1(:) ; initTheta2(:)];
nn_Thetas = initial_Thetas;


% EXECUTE COST FUNCTION

J = GXNNCostF(nn_Thetas, L1neurons, L2neurons, num_labels, X, y, lambda);


%% TRAIN NEURAL NETWORK
disp('Training Neural Network...')


options = optimset('MaxIter', maxiters);

GcostFun = @(T) GXNNCostF(T, L1neurons, L2neurons, num_labels, X, y, lambda);

[nn_Thetas, cost] = GXfmincg(GcostFun, initial_Thetas, options);



% GcostFun = @(T) GXNNCostF(T, L1neurons, L2neurons, num_labels, X, y, lambda);
% options = optimset('MaxIter', 1000,'PlotFcns',@optimplotfval);
% [nn_Thetas, cost] = fmincon(GcostFun,initial_Thetas,options);
% @optimplotx plots the current point
% @optimplotfval plots the function value
% @optimplotfunccount plots the function count (not available for fzero)




% REROLL Theta1 AND Theta2
Theta1 = reshape(nn_Thetas(1:L2neurons * (L1neurons + 1)), L2neurons, (L1neurons + 1));
Theta2 = reshape(nn_Thetas((1 + (L2neurons * (L1neurons + 1))):end), num_labels, (L2neurons + 1));

% Theta1 = reshape( nn_Thetas(1:L2neurons*L1neurons) , L2neurons, L1neurons );
% Theta2 = reshape( nn_Thetas(L2neurons*L1neurons+1:end) , num_labels, L2neurons );





%% HOW DOES NN PERFORM ON **TRAINING** SET?

[p , a , h] = GXNNpredict(Theta1, Theta2, X);

TRAINPCTCORRECT = mean(p == y);

disp('Percent accuracy on training data:')
disp( TRAINPCTCORRECT )


%% EVALUATE NEURAL NETWORK ON **TEST** SET


INPUTMX       = ADNNt(2:end,3:end);  % col 3 of ADNN is bias
INPUTYN       = ADNNt(2:end,2);

ADNNt(1:9,1:9)
INPUTMX(1:9,1:9)

INPUTMXdims   = size(INPUTMX);
INPUTYNdims   = size(INPUTYN);


% RANDOMIZE ROWS
randp         = randperm(INPUTMXdims(1));
INPUTMX       = INPUTMX(randp,:);
INPUTYN       = INPUTYN(randp,:);




X  = INPUTMX( : , UseGenes );         % INCLUDE BIAS COL

y = INPUTYN + 1;    % ADD 1 SO LABELS ARE {1,2}





[p , a , h] = GXNNpredict(Theta1, Theta2, X);

% p : prediction label {1,2}
% a : confidence level of p
% h : confidence level of p and min label

TESTPCTCORRECT = mean(p == y);

% REPORT ACCURACY FOR FULL TESTING DATASET
disp('Percent accuracy on full testing data:')
disp( TESTPCTCORRECT )

% nprtool
% load('net1.mat')
% yguess = net(X');
% mean(INPUTYN' == round(yguess))
% p = round(yguess)' + 1;

%% EXAMINE & REPORT ACCURACY FOR HIGH CONFIDENCE DATA
LowConfThresh = .51;
MidConfThresh = .80;
HiConfThresh  = .85;


numCorrect = p == y;


LowConf      = a >= LowConfThresh & a < MidConfThresh;
LowConfPCT   = (mean(p(LowConf) == y(LowConf)))*100;
LowConfNp    = numel(p(LowConf));
LowConfPCTp  = (numel(p(LowConf)) / numel(p))*100;

MidConf      = a >= MidConfThresh & a < HiConfThresh;
MidConfPCT   = (mean(p(MidConf) == y(MidConf)))*100;
MidConfNp    = numel(p(MidConf));
MidConfPCTp  = (numel(p(MidConf)) / numel(p))*100;

HiConf       = a >= HiConfThresh;                   % Logical index of high confidence predictions
HiConfPCT    = (mean(p(HiConf) == y(HiConf)))*100;  % Percent correct hi conf predictions
HiConfNp     = numel(p(HiConf));                    % Total number of hi conf predictions
HiConfPCTp   = (numel(p(HiConf)) / numel(p))*100;   % Percent of hi conf predictions
HiConfNCorr  = HiConfNp * (HiConfPCT / 100);        % Total number of correct hi conf predictions


fprintf('\nPercent accuracy on mild conf (a > %.2f) testing data: >>> %.1f%% <<<\n',...
    LowConfThresh, LowConfPCT)
fprintf('    (includes %.0f of %.0f [%.1f%%] HighConf:Total patients) \n\n',...
    LowConfNp , numel(p) , LowConfPCTp )

fprintf('\nPercent accuracy on mid conf (a > %.2f) testing data: >>> %.1f%% <<<\n',...
    MidConfThresh, MidConfPCT)
fprintf('    (includes %.0f of %.0f [%.1f%%] HighConf:Total patients) \n\n',...
    MidConfNp , numel(p) , MidConfPCTp )

fprintf('\nPercent accuracy on high conf (a > %.2f) testing data: >>> %.1f%% <<<\n',...
    HiConfThresh, HiConfPCT)
fprintf('    (includes %.0f of %.0f [%.1f%%] HighConf:Total patients) \n\n',...
    HiConfNp , numel(p) , HiConfPCTp )





%% SAVE NN RUNTIME PARAMETERS AND RESULTS



NNET.UseGenes(nn)         = {UseGenes'};
NNET.UseNgenes(nn)        = UseNgenes;
NNET.maxiters(nn)         = maxiters;
NNET.lambda(nn)           = lambda;
NNET.L2neurons(nn)        = L2neurons;
NNET.useBias(nn)          = useBias;

NNET.TRAIN_PCT(nn)        = TRAINPCTCORRECT.*100;
NNET.TEST_PCT(nn)         = TESTPCTCORRECT.*100;

NNET.TEST_LowConfPCT(nn)  = LowConfPCT;
NNET.TEST_LowConfNp(nn)   = LowConfNp;
NNET.TEST_LowConfPCTp(nn) = LowConfPCTp;

NNET.TEST_MidConfPCT(nn)  = MidConfPCT;
NNET.TEST_MidConfNp(nn)   = MidConfNp;
NNET.TEST_MidConfPCTp(nn) = MidConfPCTp;

NNET.TEST_HiConfPCT(nn)   = HiConfPCT;
NNET.TEST_HiConfNp(nn)    = HiConfNp;
NNET.TEST_HiConfPCTp(nn)  = HiConfPCTp;

NNET.UseGenes{nn}         = UseGenes;

ADNNid = ADNN(1,3:end);
UseGenesID = ADNNid(UseGenes)';
UseGenesID(1) = [];
[q,r,s] = intersect(ADGENES.VID, UseGenesID);

NNET.ADGENESrows{nn}     = r;

end
%%

% disp(' ')
% disp(NNET.UseNgenes)
% disp(NNET.maxiters)
% disp(NNET.L2neurons)
% disp(NNET.useBias)
% disp(NNET.lambda)
% disp(' ')
% disp((NNET.TEST_PCT.*100))
% disp(NNET.TEST_MidConfPCT)
% disp(NNET.TEST_HiConfPCT)
% disp((NNET.TEST_HiConfNp.*1.000001))
% disp(' ')

v0 = (NNET.TEST_HiConfNp.*(NNET.TEST_HiConfPCT ./ 100))';

v1 = [NNET.UseNgenes' NNET.maxiters' NNET.L2neurons' NNET.useBias' NNET.lambda'];

v2 = [NNET.TEST_PCT' NNET.TEST_MidConfPCT' NNET.TEST_HiConfPCT' v0 v0./5+NNET.TEST_PCT'+NNET.TEST_HiConfPCT'];



%%


RESULTS = [v1 v2];
RESULTS = [(1:size(RESULTS,1))', RESULTS];
z = sortrows(RESULTS,11);

clc; format shortG; disp(z);

% nNAN = 8;
% z(end-nNAN:end,:) = [];
% ADGENES(NNET.ADGENESrows{z(end,1)},:)


return
%% WHICH GENES WERE USED? WHICH WERE IMPORTANT?


TopRuns = flipud(z(end-99:end,1));

% ADGENES(NNET.ADGENESrows{TopRuns(1)},:)

TopRunGIDs = [];
for nn = 1:size(TopRuns,1)

    TopRunGIDs(:,nn) = ADGENES.VID(NNET.ADGENESrows{TopRuns(nn)});

end

%%
[C,ia,ic] = unique(TopRunGIDs(:));
icList = sort(ic);
[m,f] = mode(icList)
ADG = ADGENES(ADGENES.VID == C(m),:)
ADG.freq = f;

f=10;
while f>8

    icList(icList == m) = [];

    [m,f] = mode(icList);
    if f < 1; break; end

    AD = ADGENES(ADGENES.VID == C(m),:);
    AD.freq = f;

    ADG = [ADG; AD];
    disp(f)

end

ADG

%%
% [a,b,c] = intersect(ADG.VID,ADNN(1,:),'stable');
% 
% clc
% ADNN(1,c)'
% 
% [a,b,c] = intersect(xlsN,ADNN(1,:),'stable');
% ADNN(1,c)'
% 
% clc
% c
% 87615

%%
return
%% VISUALIZE NEURAL NETWORK

fprintf('\nVisualizing Neural Network... \n')

GXNNdisplayData(Theta1(:, 2:end) , 20);








%% PERFORM MANUAL GRADIENT DESCENT



[J , G] = GXNNCostFun(nn_params, input_layer_size, hidden_layer_size, num_labels, X, y, lambda);







%% TEST OTHER ADVANCED OPTIMIZATION FUNCTIONS

options = optimoptions('fminunc','Algorithm','trust-region','SpecifyObjectiveGradient',true, 'MaxIter', 100);
[optTheta, functionVal, exitFlag] = fminunc(GcostFun, initial_nn_params, options);


options = optimset('MaxIter', 100, 'PlotFcns',@optimplotfval);
[x,fval,exitflag,output] = fminsearch(GcostFun, initial_nn_params, options);



%%

num_iters = 400;

Jcost = zeros(num_iters, 1);

for iter = 1:num_iters

    theta = theta - (alpha ./ m) .* sum( (X * theta - y) .* X)';
    
    J = sum((X * theta - y).^2) ./ (2*m);

    Jcost(iter) = J;

end



% NORMAL EQUATION:
theta = pinv(X' * X) * X' * y;






% % THREE LAYER NEURAL NETWORK: random uniform {-1:1}
% L1WEIGHTS = rand( L1NEURONS , INPUTMXdims(2) ) * 2 - 1;
% L2WEIGHTS = rand( L2NEURONS , L1NEURONS ) * 2 - 1;
% L3WEIGHTS = rand( INPUTYNdims(2) , L2NEURONS  ) * 2 - 1;



% ACTIVATION FUNCTION plot(-10:.1:10,f(-10:.1:10))
% plot(  unique(INPUTMX(:),'sorted')  ,  f( unique(INPUTMX(:),'sorted') )   )

%% SET NEURAL NETWORK PARAMETERS

%-------------------------
weightScaleFactor = 0.15;
L1NEURONS = 19;
L2NEURONS = 15;
pctThresh = 93;
nEpochs = 70;
%-------------------------




%% SETUP PLOTS

% Input data from current example set
In = INPUTMX(1,:)';      % find(In)
YN = INPUTYN(1);

% In(:,2:end) = In(:,2:end) - .5;

% Propagate the signals through network
L1OUTS = f( L1WEIGHTS * In  );
L2OUTS = f( L2WEIGHTS * L1OUTS);
L3OUTS = f( L3WEIGHTS * L2OUTS);


close all; clc
fh1  = figure('Units','normalized','OuterPosition',[.2 .1 .75 .88],'Color','w');

hax1 = axes('Position',[.05 .82 .9 .14],'Color','none');
    hax1.YLim = [0 1]; hax1.YLimMode = 'manual';
    hold on;

hax2 = axes('Position',[.05 .65 .9 .14],'Color','none');
    hax2.YLim = [0 1]; hax2.YLimMode = 'manual';
    hold on;

hax3 = axes('Position',[.05 .35 .75 .25],'Color','none');
    hax3.YLim = [0 1]; hax3.YLimMode = 'manual';
    hold on;

hax4 = axes('Position',[.85 .05 .1 .5],'Color','none');
    hax4.YLim = [0 1]; hax4.YLimMode = 'manual';
    hold on;
    

hax5 = axes('Position',[.05 .05 .75 .25],'Color','none');
    hax2.YLim = [0 1]; hax2.YLimMode = 'manual';
    hold on;

axes(hax1)
ph1 = scatter( 1:length(In) , In , 500 , '.r');

axes(hax2)
ph2 = scatter( 1:length(L1OUTS) , L1OUTS , 1500 , '.r');

axes(hax3)
ph3 = scatter( 1:length(L2OUTS) , L2OUTS , 1500 , '.r');

axes(hax4)
ph4 = scatter( 1:length(L3OUTS) , L3OUTS , 1500 , '.r');


RMSmu     = zeros(nEpochs,1);

axes(hax5)
ph5 = scatter( 1:nEpochs , RMSmu , 1500 , '.r');

% t4 = text(.6,1.05,'GUESS: CASE!','Color','r','FontSize',16);    
pause(1.5)


%% -----------------------------------------------
% START MAIN TRAINING LOOP FOR 3-LAYER SUPERVISED LEARNING

% y = zeros(nEpochs,1);
figure(fh1)

nn=1; i=1;

for nn=1:nEpochs
    fprintf('\nEpoch: % .0f\n',nn)

    RMS_Err = 0;
    
    LayerB_YN = ones(INPUTMXdims(1),1);
    
    % Iterate through all examples
    for i=1:INPUTMXdims(1)
        
        % Input data from current example set
        In = INPUTMX(i,:)';
        YN = INPUTYN(i);
                
        %In(:,2:end) = In(:,2:end) - .5;
        
        % Propagate the signals through network
        L1OUTS = f(L1WEIGHTS * In    );
        L2OUTS = f(L2WEIGHTS * L1OUTS);
        L3OUTS = f(L3WEIGHTS * L2OUTS);
        
        
        
        
        delta_i = L3OUTS .* (1-L3OUTS) .* (YN-L3OUTS);           % NODE ERRORS IN LayerB
        
        delta_k = L2OUTS .* (1-L2OUTS) .* (L3WEIGHTS' * delta_i);  % NODE ERRORS IN LayerAB

        delta_j = L1OUTS .* (1-L1OUTS) .* (L2WEIGHTS' * delta_k);  % NODE ERRORS IN LayerA
                
        RMS_Err = RMS_Err + norm(YN - L3OUTS);           % ROOT MSE
        
        
        LayerB_YN(i) = L3OUTS;
                        
        % TWO LAYER NETWORK
        % U = U + weightScaleFactor .* delta_i * LayerA' ;
        % W = W + weightScaleFactor .* delta_j * In' ;
        
        % THREE LAYER NETWORK
        L2WEIGHTS = L2WEIGHTS + weightScaleFactor .* delta_k * L1OUTS' ;
        L3WEIGHTS = L3WEIGHTS + weightScaleFactor .* delta_i * L2OUTS' ;
        L1WEIGHTS = L1WEIGHTS + weightScaleFactor .* delta_j * In' ;
    end
    
    
    numCorrect = sum(INPUTYN == round(LayerB_YN));
    pctCorrect = numCorrect / INPUTYNdims(1) * 100;
    
    fprintf('\nCorrect Guesses: %.0f of %.0f  (%.2f pct.) \n', numCorrect, INPUTYNdims(1) , pctCorrect)
    
    if pctCorrect > pctThresh
        fprintf('\nPercent correct guesses is greater than: (%.2f pct.) \n', pctThresh)
        disp('Moving on to test phase...')
        pause(2)
        break
    end
    
    
    RMSmu(nn) = RMS_Err;
    
    %if i==1
    ph1.YData = In';
    ph2.YData = L1OUTS;
    ph3.YData = L2OUTS;
    ph4.YData = L3OUTS;
    ph5.YData = RMSmu;
    pause(.05)
    %end

end

% save('NN.mat')




%% RUN TEST DATA USING TRAINED NEURAL NETWORK ABOVE
clearvars -except L1WEIGHTS L2WEIGHTS L3WEIGHTS ADNN ADNNt ADMx ADGENES
close all;

disp(' ------------------------ ')
disp('      TEST SESSION        ')
disp(' ------------------------ ')


ADNNt(2:end,3) =  1; % BIAS COL#3:
INPUTMX        = ADNNt(2:end,3:end);
INPUTYN        = ADNNt(2:end,2);
% SNPSITE        = SNPSITE - .5;

INPUTMXdims    = size(INPUTMX);
INPUTYNdims    = size(INPUTYN);

randp         = randperm(INPUTMXdims(1));
INPUTMX       = INPUTMX(randp,:);
INPUTYN       = INPUTYN(randp,:);

disp([INPUTYN(1:6,:) INPUTMX(1:6,1:8)])


f = @(x) (1./(1+exp(-x)));


for nn=1:1

    RMS_Err = 0;
    
    LayerB_YN = ones(INPUTMXdims(1),1);
    
    for i=1:INPUTMXdims(1)
        
        In = INPUTMX(i,:)';
        YN = INPUTYN(i);
        
        L1OUTS = f(L1WEIGHTS*In);
        L2OUTS = f(L2WEIGHTS*L1OUTS);
        L3OUTS = f(L3WEIGHTS*L2OUTS);
        delta_i = L3OUTS .* (1-L3OUTS) .* (YN-L3OUTS);
        delta_k = L2OUTS .* (1-L2OUTS) .* (L3WEIGHTS' * delta_i);
        delta_j = L1OUTS .* (1-L1OUTS) .* (L2WEIGHTS' * delta_k);
        RMS_Err = RMS_Err + norm(YN - L3OUTS);
        LayerB_YN(i) = L3OUTS;
        
    end
    
    numCorrect = sum(INPUTYN == round(LayerB_YN));
    pctCorrect = numCorrect / INPUTYNdims(1) * 100;
    
    fprintf('Correct Guesses: %.0f of %.0f  (%.1f%%) \n', ...
            numCorrect, INPUTYNdims(1) , pctCorrect)
    
end

fprintf('Number of CASES: % .0f\n',sum(INPUTYN))
fprintf('Number of CTRLS: % .0f\n',INPUTYNdims(1) - sum(INPUTYN))





%% TEST IDENTITY MATRIX TO DETERMINE WHAT SNPS IMPACT NN OUTPUT MOST

clearvars -except L1WEIGHTS L2WEIGHTS L3WEIGHTS ADNN ADNNt ADMx ADGENES


INPUTMX       = ADNN(2:end,3:end);
INPUTYN       = ADNN(2:end,2);

INPUTMXdims   = size(INPUTMX);
INPUTYNdims   = size(INPUTYN);


IdentityMx = eye(INPUTMXdims(2));
IdentityMx(:,1) = 1;


OutputLayer_Response = zeros(INPUTMXdims(2),1);

f = @(x) (1./(1+exp(-x)));
 
for i=1:INPUTMXdims(2)
    In = IdentityMx(i,:)';
    L1OUTS = f(L1WEIGHTS*In);
    L2OUTS = f(L2WEIGHTS*L1OUTS);
    L3OUTS = f(L3WEIGHTS*L2OUTS);
    OutputLayer_Response(i) = L3OUTS;
end

% [NNcoRes,NNcoInd] = sort(OutputLayer_Response(2:end));
[NNcoRes,NNcoInd] = sort(OutputLayer_Response);
NNcaRes = flipud(NNcoRes);
NNcaInd = flipud(NNcoInd);

format shortG
[NNcaRes(1:10) NNcaInd(1:10) NNcoRes(1:10) NNcoInd(1:10)]

Ng = 20;

CTRLgenes = NNcoInd(1:Ng);
CTRLgeneWeights = NNcoRes(1:Ng);

CASEgenes = NNcaInd(1:Ng);
CASEgeneWeights = NNcaRes(1:Ng);




% INPUTMX = ADNN(2:end,3:end);
INPUTMX = ADNN(:,3:end);


CASEVID = INPUTMX(1, NNcaInd);
CTRLVID = INPUTMX(1, NNcoInd);

CASEWt = NNcaRes;
CTRLWt = NNcoRes;

load('GXUNDATA.mat','UNMX')


UNMX( CASEVID(1:50) , :)
UNMX( CTRLVID(1:50) , :)


[CAv, CAi, CAj] = intersect(ADGENES.VID, CASEVID(1:50)' );
[COv, COi, COj] = intersect(ADGENES.VID, CTRLVID(1:50)' );

CAGENES = ADGENES(CAi,:); 
CAGENES.CACO = ones(size(CAGENES,1),1);

COGENES = ADGENES(COi,:); 
COGENES.CACO = zeros(size(COGENES,1),1);

CAGENES.Wt = CASEWt(CAj);
COGENES.Wt = CTRLWt(COj);

CACOGENES = [CAGENES; COGENES];



keepvars = {...
'keepvars'
'CASEVID'
'CTRLVID'
'CASEVIDWt'
'CTRLVIDWt'
'ADMX'
'L1WEIGHTS'
'L2WEIGHTS'
'L3WEIGHTS'
'ADNN'
'ADNNt' 
'ADMx'
'OutputLayer_Response'
'CASEgeneWeights'
'CTRLgeneWeights'
'CASEgenes'
'CTRLgenes'
'ADGENES'
};
clearvars('-except',keepvars{:})

%%



INPUTX = ADNN(1,3:end);

% IF BIAS IS STILL IN OUTPUT LAYER SOMETHING IS WRONG
% CASEgeneWeights(CASEgenes == 1) = [];
% CTRLgeneWeights(CTRLgenes == 1) = [];
% CASEgenes(CASEgenes == 1) = [];
% CTRLgenes(CTRLgenes == 1) = [];


[CAv, CAi, CAj] = intersect(ADGENES.VID, INPUTX(1,CASEgenes'));
[COv, COi, COj] = intersect(ADGENES.VID, INPUTX(1,CTRLgenes'));

ADGENES(CAi,:)
ADGENES(COi,:)



load('GXUNDATA.mat','UNMX')

UNMX(CASEgenes , :)
UNMX(CTRLgenes , :)


%%

% GRAPH AND TABLE
% lintrans = @(x,a,b,c,d) (c.*(1-(x-a)./(b-a)) + d.*((x-a)./(b-a)));

NNresp = OutputLayer_Response(2:end);

close all
fh1  = figure('Units','normalized','OuterPosition',[.02 .06 .95 .9],'Color','w','MenuBar','none');
hax1 = axes('Position',[.05 .05 .90 .90],'Color','none','XTick',[],'YTick',[]);

plot(NNresp)
hold on
% scatter(1:numel(NNresp),NNresp,200,'.b')
% hold on
scatter(CASEgenes, CASEgeneWeights, 900,'.r')
hold on
scatter(CTRLgenes, CTRLgeneWeights, 900,'.g')



return


%% RUN PRINCIPLE COMPONENTS ANALYSIS
%{
Ncomps = 1000;

opt = statset('pca');
opt.MaxIter = 5000;

[PCAScof,PCASval,~,~,PCASexp,~] = pca(...
    NNOBSca,'Options',opt,'Algorithm','svd','NumComponents',Ncomps,'Centered',false);

[PCONcof,PCONval,~,~,PCONexp,~] = pca(...
    NNOBSco,'Options',opt,'Algorithm','svd','NumComponents',Ncomps,'Centered',false);

% PCAStopcof = PCAScof(1,:);
% PCAScentered = PCASval*PCAScof';
% PCAStsredu = mahal(PCASval,PCASval);
% PCAStsqdiscard = PCASts - PCAStsredu;


clc
sum(PCASexp(1:100))
sum(PCONexp(1:100))

% format shortG
% PCASval(1:5 , :)


%% PREP PCA DATA BACK INTO NEURAL NETWORK SPARSE MATRIX

PCACO = [PCASval'; PCONval'];

PCANN = padarray(PCACO,[1 3],0,'pre');

nrows = size(PCANN,1);
ncols = size(PCANN,1);


PCANN(2:Ncomps+1 , 2)  =  1;                % COL#2: case/ctrl 1/0
PCANN(2:end , 3)       =  1;                % COL#3: bias
% PCANN(1, 4:end)      =  VID';             % ROW#1: VID
% PCANN(2:end, 1)      =  [CaseID; CtrlID]; % COL#1: subject IDs

%}

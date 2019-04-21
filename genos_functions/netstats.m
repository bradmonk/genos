function [varargout] = netstats(net,TRAINMX,TRAINLAB,TESTMX,TESTLAB,Pfilter,HICI,LOCI,doROC)


clc;


TRvec = net(TRAINMX);
TEvec = net(TESTMX);

% vec2ind can be used as a shortcut to get row index of max value in each column
TRclass = vec2ind(TRvec);
TEclass = vec2ind(TEvec);


TRGUESS = round(TRvec);
TEGUESS = round(TEvec);

TRCASEn = sum(TRAINLAB(1,:));
TRCTRLn = sum(TRAINLAB(2,:));
TECASEn = sum(TESTLAB(1,:));
TECTRLn = sum(TESTLAB(2,:));

disp(' ')
fprintf('%4.0f  number of training cases\n',TRCASEn)
fprintf('%4.0f  number of training ctrls\n',TRCTRLn)
fprintf('%4.0f  number of testing cases\n',TECASEn)
fprintf('%4.0f  number of testing ctrls\n',TECTRLn)
fprintf('%4.0f  number of times net guessed case on training data\n',sum(TRGUESS(1,:)))
fprintf('%4.0f  number of times net guessed ctrl on training data\n',sum(TRGUESS(2,:)))
fprintf('%4.0f  number of times net guessed case on testing data\n',sum(TEGUESS(1,:)))
fprintf('%4.0f  number of times net guessed ctrl on testing data\n',sum(TEGUESS(2,:)))
disp(' ')



yguess = net(TRAINMX);
HITRPCTALL = mean(mean(TRAINLAB == round(yguess))) * 100;
disp(' '); disp('------------ TRAINING DATA ------------')
fprintf('%.2f  Percent correct on all training data \n',HITRPCTALL)


yi = yguess>HICI | yguess<LOCI;
yi = yi(1,:);
HITRPCT = mean(mean(  TRAINLAB(:,yi) == round(yguess(:,yi))  )) * 100;
HITRCON = (sum(yi) / numel(yi)) * 100;
fprintf('%.2f  Percent correct on high-confidence training data \n',HITRPCT)
fprintf('%.2f  Percent of training data registered high-confidence \n',HITRCON)
disp(' ');


yguess = net(TESTMX);
HITEPCTALL = mean(mean(TESTLAB == round(yguess))) * 100;
disp(' '); disp('------------- TEST DATA ---------------')
PERF.all = sprintf('%.2f  Percent correct on all hold-out test data',HITEPCTALL);
disp(PERF.all)









HIC_TOPTHRESH = HICI;
HIC_BOTTHRESH = LOCI;

MIDC_TOPTHRESH = .70;
MIDC_BOTTHRESH = .30;





yi = yguess>MIDC_TOPTHRESH | yguess<MIDC_BOTTHRESH;
yi = yi(1,:);

MIDTEPCT = mean(mean(  TESTLAB(:,yi) == round(yguess(:,yi))  )) * 100;
MIDTEPOP = (sum(yi) / numel(yi)) * 100;

disp(' '); disp('--- MEDIUM CONFIDENCE ---')
PERF.midpct = sprintf('%.2f  Percent correct on mid-confidence hold-out test data',MIDTEPCT);
PERF.midpop = sprintf('%.2f  Percent of hold-out test data registered mid-confidence',MIDTEPOP);
disp(PERF.midpct);disp(PERF.midpop)











yi = yguess>HIC_TOPTHRESH | yguess<HIC_BOTTHRESH;
yi = yi(1,:);

HITEPCT = mean(mean(  TESTLAB(:,yi) == round(yguess(:,yi))  )) * 100;
HITEPOP = (sum(yi) / numel(yi)) * 100;

disp(' '); disp('--- HIGH CONFIDENCE ---')
PERF.hipct = sprintf('%.2f  Percent correct on high-confidence hold-out test data',HITEPCT);
PERF.hipop = sprintf('%.2f  Percent of hold-out test data registered high-confidence',HITEPOP);
disp(PERF.hipct);disp(PERF.hipop)

disp(' ');
fprintf('**high-confidence definition: %.2f > activation > %.2f \n',HIC_BOTTHRESH,HIC_TOPTHRESH)
fprintf('**medium-confidence definition: %.2f > activation > %.2f \n',MIDC_BOTTHRESH,MIDC_TOPTHRESH)
disp(' ')












if doROC > 0
%------------------------------------
% PLOT HISTOGRAMS AND ROC CURVE

close all;
fh1=figure('Units','normalized','OuterPosition',[.10 .06 .7 .9],'Color','w','MenuBar','none');
h1 = axes('Position',[.06 .62 .40 .33],'Color','none');
h2 = axes('Position',[.53 .62 .40 .33],'Color','none','XDir','reverse'); box on; hold on;
h3 = axes('Position',[.06 .06 .40 .40],'Color','none'); axis off; hold on;
h4 = axes('Position',[.06 .06 .40 .40],'Color','none','XDir','reverse'); axis off; hold on;
h5 = axes('Position',[.53 .06 .40 .45],'Color','none'); hold on;
h6 = axes('Position',[.06 .46 .40 .10],'Color','none'); axis off; hold on;


axes(h1); 
p1=histogram(yguess(1,TESTLAB(1,:)==1),20,'BinLimits',[0 1],...
        'EdgeColor','k','FaceColor','r'); title('CASE CONFIDENCE'); box on;

axes(h2); 
p2=histogram(yguess(2,TESTLAB(1,:)==0),20,'BinLimits',[0 1]); 
        title('CTRL CONFIDENCE');  box on;

axes(h3); 
p3=histogram(yguess(1,TESTLAB(1,:)==1),20,'BinLimits',[0 1], ...
        'EdgeColor','k','FaceColor','r');  box off;

axes(h4); 
p4=histogram(yguess(2,TESTLAB(1,:)==0),20,'BinLimits',[0 1]); 
        xlabel('\leftarrow CTRLS vs CASE \rightarrow');


% Add lines to histogram plots
ymax = max([h1.YLim,h2.YLim,h3.YLim,h4.YLim]);
h1.YLim = [0 ymax];
h2.YLim = [0 ymax];
h3.YLim = [0 ymax];
h4.YLim = [0 ymax];

axes(h1); line([HICI HICI],[0 ymax],'Color','k','LineWidth',1,'LineStyle','--')
axes(h1); line([LOCI LOCI],[0 ymax],'Color','k','LineWidth',1,'LineStyle','--')
axes(h2); line([HICI HICI],[0 ymax],'Color','k','LineWidth',1,'LineStyle','--')
axes(h2); line([LOCI LOCI],[0 ymax],'Color','k','LineWidth',1,'LineStyle','--')
axes(h4); line([.5 .5],[0 ymax],'Color','k','LineWidth',1,'LineStyle','--')


% Add summary stats to figure
axes(h6); plot(1,1)
PERF.pvn = ['P(' num2str(Pfilter) ')=' num2str(size(TRAINMX,1)) ' variants'];
text(0,2.0, PERF.pvn,'FontSize',12); 
text(0,1.5, PERF.all,'FontSize',12); 
text(0,1.0, PERF.hipct,'FontSize',12); 
text(0,0.5, PERF.hipop,'FontSize',12)

%PERF.mix = sprintf('%.0f%% scrambled training data ',MIXRATIO*100);
%text(1.5,2.0, PERF.mix,'FontSize',12)




% Add colored patch overlay to combined histogram
axes(h3);
x1 = [LOCI LOCI HICI HICI];
y1 = [0 ymax/1 ymax/1 0];
x2 = [0 0 LOCI LOCI];
y2 = [0 ymax/1 ymax/1 0];
x3 = [HICI HICI 1 1];
y3 = [0 ymax/1 ymax/1 0];
patch(x1,y1,'black','FaceAlpha',.1);
patch(x2,y2,'green','FaceAlpha',.1);
patch(x3,y3,'red','FaceAlpha',.1);


% PLOT ROC CURVE
simTRAIN = sim(net,TRAINMX);
simTEST = sim(net,TESTMX);



[TPR,FPR,THR] = roc(TRAINLAB,simTRAIN);
if size(TPR{1},2) == size(TPR{2},2)
    TPR = cell2mat(TPR');
    FPR = cell2mat(FPR');
    THR = cell2mat(THR');

    axes(h5)
    ph1 = plot(FPR',TPR','LineWidth',5);
    ph2 = line([0 1],[0 1],'Color','k','LineStyle','--','LineWidth',1);
    box on; grid on

else
    TP = TPR{1}';
    FP = FPR{1}';
    TH = THR{1}';

    axes(h5)
    ph1 = plot(FP',TP','LineWidth',5);
    ph2 = line([0 1],[0 1],'Color','k','LineStyle','--','LineWidth',1);
    box on; grid on; hold on

    TP = TPR{2}';
    FP = FPR{2}';
    TH = THR{2}';

    plot(FP',TP','LineWidth',5);
    line([0 1],[0 1],'Color','k','LineStyle','--','LineWidth',1);
    box on; grid on
end


% END IF doROC
end







NNPERF = [HITRPCTALL,HITEPCTALL,...
          MIDTEPCT,MIDTEPOP,...
          HITEPCT,HITEPOP];



if nargout==1
    varargout = {NNPERF};
end

end


%% DEEP LEARNING: CONVOLUTIONAL NEURAL NETWORK
%{
clc
clearvars -except ADMX SNPCASE SNPCTRL PHE CASEID CTRLID CASEMX CTRLMX...
COHORTS AMX AMXCASE AMXCTRL ADNN caMX coMX Pfilter...
TRAINX TRAINL TESTX TESTL TRAINMX TRAINLAB TESTMX TESTLAB net

load lettersTrainSet

TRAINMX  = TRAINX';
TRAINLAB = TRAINL';

TESTMX  = TESTX';
TESTLAB = TESTL';


[TTACTUAL,~,~] = find(TRAINLAB);
zTRAINLAB = categorical(TTACTUAL);


zTRAINMX = TRAINMX;
sz = size(zTRAINMX,1);
re = sz-1800;
rp = randperm(sz,re);
zTRAINMX(rp,:) = [];
zTRAINMX = reshape(zTRAINMX,50,36,1,size(zTRAINMX,2));



clc
size(zTRAINMX)
size(zTRAINLAB)
size(XTrain)
size(TTrain)


layers = [imageInputLayer([50 36 1]);
          convolution2dLayer(3,3);
          reluLayer();
          maxPooling2dLayer(2,'Stride',2);
          fullyConnectedLayer(size(TRAINLAB,1));
          softmaxLayer();
          classificationLayer()];

% convolution2dLayer(5,16);
% options = trainingOptions('sgdm');

options = trainingOptions('sgdm',...
    'LearnRateSchedule','piecewise',...
    'LearnRateDropFactor',0.2,...
    'LearnRateDropPeriod',5,...
    'MaxEpochs',20,...
    'MiniBatchSize',64,...
    'Verbose',false,...
    'Plots','training-progress',...
    'OutputFcn',@(info)stopIfAccuracyNotImp(info,3));


net = trainNetwork(zTRAINMX,zTRAINLAB,layers,options);



[TTACTUAL,~,~] = find(TESTLAB);
zTESTLAB = categorical(TTACTUAL);


zTESTMX = TESTMX;
sz = size(zTESTMX,1);
re = sz-1800;
rp = randperm(sz,re);
zTESTMX(rp,:) = [];
zTESTMX = reshape(zTESTMX,50,36,1,size(zTESTMX,2));

YTest = classify(net,zTESTMX);

accuracy = sum(YTest == zTESTLAB)/numel(zTESTLAB);

disp(accuracy)



%% DEEP LEARNING: LSTM NETWORK
%}
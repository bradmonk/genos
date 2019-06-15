function [CONMX, MU, varargout] = svcmx(PVMX, LABS, ACT, SVM, XLOCI)
%======================================================================
% BUILD CONFUSION MATRICES
%======================================================================
%  %CORRECT     NCASE_ACTUAL   NCTRL_ACTUAL   %CASE_ACTUAL
%  NCASE_PRED       TP             FP             PPV
%  NCTRL_PRED       FN             FN             NPV
%  %CASE_PRED      TPR            TNR              J
%------------------------------------------------------------------


    [ERR, CMX,~, CMR] = confusion(LABS,ACT);
    COR = 1-ERR;
    CMX = CMX';

    


    CONMX = nan(4,4);

    CONMX(2:3,2:3) = CMX;

    CONMX(2:3,1) = [CMX(1,1)+CMX(1,2);
                    CMX(2,1)+CMX(2,2)];

    CONMX(1,2:3) = [CMX(1,1)+CMX(2,1),...
                    CMX(1,2)+CMX(2,2)];


    % CHECK IF PREDICTION COUNTS MATCH
    %all(sum(round(ACT),2) == CONMX(2:3,1))

    % CHECK IF ACTUAL LABEL COUNTS MATCH
    %all(sum(LABS,2)' == CONMX(1,2:3))


    TPR = CONMX(2,2) / CONMX(1,2);
    TNR = CONMX(3,3) / CONMX(1,3);
    PPV = CONMX(2,2) / CONMX(2,1);
    NPV = CONMX(3,3) / CONMX(3,1);

    PACASE = CONMX(1,2) / (CONMX(1,2) + CONMX(1,3));
    PPCASE = CONMX(2,1) / (CONMX(2,1) + CONMX(3,1));



%======================================================================
% BUILD CONFUSION MATRICES
%======================================================================
%  %CORRECT     NCASE_ACTUAL   NCTRL_ACTUAL   %CASE_ACTUAL
%  NCASE_PRED       TP             FP             PPV
%  NCTRL_PRED       FN             FN             NPV
%  %CASE_PRED      TPR            TNR              J
%------------------------------------------------------------------
    CONMX(1,4) = PACASE;
    CONMX(4,1) = PPCASE;
    CONMX(4,2) = TPR;
    CONMX(4,3) = TNR;
    CONMX(2,4) = PPV;
    CONMX(3,4) = NPV;
    CONMX(1,1) = COR;
    CONMX(4,4) = TPR + TNR - 1;

    



%------------------------------------------------------------------



    HACT = ACT - .5;


    % ALL: SPLIT PREDICTIONS FOR CASES AND CONTROLS
    ACT_CASE = HACT(1,LABS(1,:)==1) .*  1;
    ACT_CTRL = HACT(1,LABS(1,:)==0) .* -1;


    % ALL: DESCRETIZE PREDICTIONS
    TF_CASE = ACT_CASE >  0;
    TF_CTRL = ACT_CTRL >= 0;


    % ALL: PROPORTION CORRECT
    MU_CASE = nanmean(TF_CASE);
    MU_CTRL = nanmean(TF_CTRL);

    

    % ALL/HIGH: HAS HIGH ACTIVATION?
    HIQ_CASE = quantile(abs(ACT_CASE),.75);
    HIQ_CTRL = quantile(abs(ACT_CTRL),.75);

    ISHI_CASE = abs(ACT_CASE) > HIQ_CASE;
    ISHI_CTRL = abs(ACT_CTRL) > HIQ_CTRL;


    % HIGH: PROPORTION CORRECT
    HIMU_CASE = nanmean(ACT_CASE(ISHI_CASE) > 0);
    HIMU_CTRL = nanmean(ACT_CTRL(ISHI_CTRL) > 0);


    MU = [MU_CASE MU_CTRL HIMU_CASE HIMU_CTRL];





%======================================================================
%% BUILD CONFUSION MATRICES
%======================================================================
keyboard
%%





%--------------------------------------------------------------------------
% TRAIN SUPPORT VECTOR MACHINE
%--------------------------------------------------------------------------
% ML_predictions(SVMOD,PVTR(:,10:end),PVTR(:,2))


% svc = fitcsvm(PVTR(:,10:end),PVTR(:,2),'ClipAlphas',true);




% imp = SVM.predictorImportance;

SupportVectors      = SVM.SupportVectors;
SupportVectorLabels = SVM.SupportVectorLabels;


MX = SupportVectors .* SupportVectorLabels;
IMPORTANCE = sum(MX);


[IMPORTANCE , ix] = sort(IMPORTANCE,'descend')
GENENAME = XLOCI.GENE(ix);



%LOOPDATA.SVM(:,1:2,IJ) = IMPORTANCE;





% DETERMINE & PLOT IMPORTANT SVM FEATURES
%--------------------------------------------------------------------------
clc; close all;
fh1 = figure('Units','normalized','Position',[.1 .1 .55 .7],'Color','w');
ax1 = axes('Position',[.07 .15 .90 .8],'Color','none');
ph1 = bar(IMPORTANCE);
ax1.XTick = 1:numel(IMPORTANCE);
ax1.XTickLabel = GENENAME;
ax1.XTickLabelRotation = 90;
ax1.FontSize = 14;
ylim([-40 40]);
title('Feature importance');
ylabel('Support Vectors');
% pause(.5); 
% print(['/Users/bradleymonk/Desktop/APOE_34_44_' num2str(IJ)],'-dpng','-r0'); 
% pause(.5);











%--------------------------------------------------------------------------
%% JUNK TEST



[ranks,weights] = relieff(PVMX(:,10:end),PVMX(:,2),10);
bar(weights(ranks))
xlabel('Predictor rank')
ylabel('Predictor importance weight')




[ranks,weights] = relieff(PVMX(:,10:end),PVMX(:,2),10,'categoricalx','on')
bar(weights(ranks))
xlabel('Predictor rank')
ylabel('Predictor importance weight')







clc;
mdl = fscnca( PVMX(:,10:end) , PVMX(:,2) ,'Verbose',1);



cvp = cvpartition(PVMX(:,2),'holdout',200)
Xtrain = PVMX(cvp.training,10:end);
ytrain = PVMX(cvp.training,2);
Xtest  = PVMX(cvp.test,10:end);
ytest  = PVMX(cvp.test,2);

nca = fscnca(Xtrain,ytrain,'FitMethod','none');
L = loss(nca,Xtest,ytest)
nca = fscnca(Xtrain,ytrain,'FitMethod','exact','Lambda',0,...
      'Solver','sgd','Standardize',true);
L = loss(nca,Xtest,ytest)
cvp = cvpartition(ytrain,'kfold',5);
numvalidsets = cvp.NumTestSets;
n = length(ytrain);
lambdavals = linspace(0,20,20)/n;
lossvals = zeros(length(lambdavals),numvalidsets);
for i = 1:length(lambdavals)
    for k = 1:numvalidsets
        X = Xtrain(cvp.training(k),:);
        y = ytrain(cvp.training(k),:);
        Xvalid = Xtrain(cvp.test(k),:);
        yvalid = ytrain(cvp.test(k),:);

        nca = fscnca(X,y,'FitMethod','exact', ...
             'Solver','sgd','Lambda',lambdavals(i), ...
             'IterationLimit',30,'GradientTolerance',1e-4, ...
             'Standardize',true);
                  
        lossvals(i,k) = loss(nca,Xvalid,yvalid,'LossFunction','classiferror');
    end
end

close all; figure
ph=plot(mdl.FeatureWeights,'.','MarkerSize',30)
grid on
xlabel('Feature index')
ylabel('Feature weight')
meanloss = mean(lossvals,2);
figure()
plot(lambdavals,meanloss,'ro-')
xlabel('Lambda')
ylabel('Loss (MSE)')
grid on

[~,idx] = min(meanloss) % Find the index
bestlambda = lambdavals(idx) % Find the best lambda value
bestloss = meanloss(idx)

nca = fscnca(Xtrain,ytrain,'FitMethod','exact','Solver','sgd',...
    'Lambda',bestlambda,'Standardize',true,'Verbose',1);


close all; figure()
plot(nca.FeatureWeights,'ro')
xlabel('Feature index')
ylabel('Feature weight')
grid on


tol    = 0.02;
selidx = find(nca.FeatureWeights > tol*max(1,max(nca.FeatureWeights)))

L = loss(nca,Xtest,ytest)


features = Xtrain(:,selidx);

svmMdl = fitcsvm(features,ytrain);

L = loss(svmMdl,Xtest(:,selidx),ytest)









%%
[PCAc,PCAv,a,b,c,d] = pca(PVMX(:,10:end));


X = PCAv(:,1:20);
Y = PVMX(:,2);




SIGFUN = @(g,u,v,c) tanh(g*u*v' + c);

SVMOD = fitcsvm(X,Y,'KernelFunction','@SIGFUN','Standardize',true);


clc; close all; figure
gscatter(  X(:,1),  X(:,2), Y ,[],[],20);
legend('Location','Northwest'); axis tight






[x1Grid,x2Grid] = meshgrid(min(X(:,1)):.02:max(X(:,1)), min(X(:,2)):.02:max(X(:,2)));

xGrid = [x1Grid(:),x2Grid(:)];

[~,scores1] = predict(SVMOD,xGrid);


close all; figure;
h(1:2) = gscatter(X(:,1),X(:,2),Y);
hold on
h(3) = plot(X(Mdl1.IsSupportVector,1), X(Mdl1.IsSupportVector,2),'ko','MarkerSize',10);

    
% Support vectors
contour(x1Grid,x2Grid,reshape(scores1(:,2),size(x1Grid)),[0 0],'k');

title('Scatter Diagram with the Decision Boundary')
legend({'-1','1','Support Vectors'},'Location','Best'); hold off













%% ------------------------------------------------------------------------






% close all
% h1 = gscatter(...
%     (PVTR(:,2)-.5) .* rand(size(PVTR(:,2))),...
%     PVTR(:,10:end) .* rand(size(PVTR(:,10:end))),...
%     PVTR(:,5),'rb','v^',[],'off');





%%



SVMModel = fitcsvm(TR_MX,TR_VL,'ClassNames',{'CTRL','CASE'},'Standardize',true);


ScoreSVMModel = fitPosterior(svm);
[label,scores] = resubPredict(svm);
[~,postProbs] = resubPredict(ScoreSVMModel);
table(TR_VL(1:10),label(1:10),scores(1:10,2),postProbs(1:10,2),'VariableNames',...
{'TrueLabel','PredictedLabel','Score','PosteriorProbability'})

mean(TR_VL==label)
mean(TR_VL==label)



SVMModel = fitcsvm(X,Y,'ClassNames',{'b','g'},'Standardize',true);




%%
close all; hold on
h1 = gscatter(SL,SW,group,'rb','v^',[],'off');
hold on
[X,Y] = meshgrid(linspace(4.5,8),linspace(2,4));
X = X(:); 
Y = Y(:);
[C,err,P,logp,coeff] = classify([X Y],[SL SW],group,'Quadratic');
gscatter(X,Y,C,'rb','.',1,'off');
K = coeff(1,2).const;
L = coeff(1,2).linear; 
Q = coeff(1,2).quadratic;
f = @(x,y) K + L(1)*x + L(2)*y + Q(1,1)*x.*x + (Q(1,2)+Q(2,1))*x.*y + Q(2,2)*y.*y;
h2 = fimplicit(f,[4.5 8 2 4]);
set(h2,'Color','m','LineWidth',2,'DisplayName','Decision Boundary')
axis([4.5 8 2 4])
xlabel('Sepal Length')
ylabel('Sepal Width')











%--------------------------------------------------------------------------
%%          HISTOGRAM AREA BAR PLOTS
%--------------------------------------------------------------------------

    close all;
    fh1 = figure('Units','pixels','Position',[10 35 760 760],'Color','w');
    edges = (-.5:.02:.5);
    ax1=subplot(3,7,1:3); histogram(ACT,(0:.02:1));     title('Neural Net Activation')
    ax2=subplot(3,7,5:7); histogram(HACT,edges);        title('Classifier Zeroed Cutoff')
    ax3=subplot(3,7,8:10); histogram(ACT_CASE,edges);    title('Case Activation')
    ax4=subplot(3,7,12:14); histogram(ACT_CTRL.*-1,edges);title('Control Activation')

    ax5=subplot(3,7,20:21); bar(([HIMU_CASE,HIMU_CTRL].*100),.5,'FaceColor',[.2 .2 .5]);
        hold on; bar(([MU_CASE,MU_CTRL].*100),.25,'FaceColor',[.9 .1 .1]);
    grid on; ylabel('Pct. Correct')
    legend({'Top Qtr.','All Qtr.'},'Location','Northwest','NumColumns',2)
    ax5.YLim = [0 116]; 
    ax5.YTick = [0 25 50 75 100]; 
    ax5.XTickLabels = {'Case','Control'};
    title('Performance Summary')

    ax6=subplot(3,7,15:18); 
        h1=histogram((ACT_CASE(ACT_CASE>0)),edges); hold on;
        h2=histogram((ACT_CTRL(ACT_CTRL>0).*-1),edges); hold on;
        h3=histogram((ACT_CASE(ACT_CASE<0)),edges); hold on;
        h4=histogram((ACT_CTRL(ACT_CTRL<0).*-1),edges); hold on;
        title('CTRL Activation * -1')
        legend({'Case Correct','Ctrl Correct',...
                'Case Incorrect','Ctrl Incorrect'},...
                'Location','best','NumColumns',2);
        title('Classifier Performance Histogram');


    hold off
    x = mean([edges(1:end-1); edges(2:end)]);
    ax6=subplot(3,7,15:18); 
        XAREA = x';
        YAREA = [h3.Values; h4.Values;h1.Values; h2.Values]';
        area(XAREA,YAREA)
        legend({'Case Miss','Ctrl Miss',...
                'Case Hit','Ctrl Hit'},...
                'Location','best');
        title('Classifier Performance Area');


    XYAREA = [XAREA YAREA];
    varargout = {XYAREA};

%--------------------------------------------------------------------------
%%


end
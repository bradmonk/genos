function [varargout] = regressionresults(ML,TRX,TEX,LOCI)
%% predictionresults

% keyboard


TRAINX = TRX(:,10:end);
TRAINL = TRX(:,2);

TESTX = TEX(:,10:end);
TESTL = TEX(:,2);


% Yreal = TRX(:,9);
% Ypred = ML.predictFcn(TRAINX);
% Labs = TRAINL;



Yreal = TEX(:,9);
Ypred = ML.predictFcn(TESTX);
Labs = TESTL;



%################   TWO PACK   ################
close all
fh1 = figure('Units','normalized','OuterPosition',[.01 .06 .5 .7],...
    'Color','w','MenuBar','none');
ax1 = axes('Position',[.08 .09 .82 .83],'Color','none'); hold on;
ax2 = axes('Position',[.08 .09 .82 .83],'Color','none'); hold on;

axes(ax1); colormap(ax1,[.5 .5 .5; .5 .5 .5])
ph1=scatter(1:numel(Yreal),Yreal,30,Labs,'filled');
ph1.MarkerFaceAlpha = .8;

axes(ax2); colormap(ax2,[.1 .99 .3; .99 .4 .4])
ph2=scatter(1:numel(Ypred),Ypred,50,Labs,'filled');
ph2.MarkerFaceAlpha = .7;

YL = [ax1.YLim; ax2.YLim];
ax1.YLim = [min(YL(:)) max(YL(:))]; %ax1.YLim = [-170,170];
ax2.YLim = [min(YL(:)) max(YL(:))]; %ax2.YLim = [-170,170];


figure
boxplot(Ypred,Labs,'PlotStyle','compact')


FIT = ML.predictFcn(TESTX);
GUESSES = ((FIT < 0) .*2 -1) == TESTL;
PCTCORRECT = mean(GUESSES) .* 100;
[LABgi,~] = findgroups(TESTL);
[FITgi,~] = findgroups(GUESSES);
fields = fieldnames(ML);
MLname = fields{2};






disp(' ');
disp('--------- HOLD-OUT TEST SET RESULTS ------------')
disp(MLname);
fprintf('%4.0f   number of test  case\n',sum(LABgi==2))
fprintf('%4.0f   number of test  ctrl\n',sum(LABgi==1))
fprintf('%4.0f   times predicted case\n',sum(FITgi==2))
fprintf('%4.0f   times predicted ctrl\n',sum(FITgi==1))
fprintf('%5.1f  percent correct\n',PCTCORRECT)
disp(' ')



FIT = ML.predictFcn(TESTX);
F = (FIT > 40) | (FIT < -40);

GUESSES = ((FIT(F) < 0) .*2 -1) == TESTL(F);
PCTCORRECT = mean(GUESSES) .* 100;
[LABgi,~] = findgroups(TESTL(F));
[FITgi,~] = findgroups(GUESSES);
fields = fieldnames(ML);
MLname = fields{2};

disp(' ');
disp('--------- HOLD-OUT TEST SET RESULTS ------------')
disp(MLname);
fprintf('%4.0f   number of test  case\n',sum(LABgi==2))
fprintf('%4.0f   number of test  ctrl\n',sum(LABgi==1))
fprintf('%4.0f   times predicted case\n',sum(FITgi==2))
fprintf('%4.0f   times predicted ctrl\n',sum(FITgi==1))
fprintf('%5.1f  percent correct\n',PCTCORRECT)
disp(' ')


% plotPartialDependence(ML.RegressionSVM,1,'Conditional','centered')
% plotPartialDependence(ML.RegressionSVM,3,'Conditional','absolute')
% pt1 = linspace(-2,6,20)';
% pt2 = linspace(-2,6,20)';
% ax = plotPartialDependence(ML.RegressionSVM,[1,3],'QueryPoints',[pt1 pt2]);
% view(140,30) % Modify the viewing angle


%%
keyboard
%% ################   PREDICTOR IMPORTANCE   ################

% ng = numel(LOCI.GENE);
% ordr   = string((1:ng))';
% GHEAD = strcat( repmat("i",ng,1), ordr  , repmat("_",ng,1) , LOCI.GENE );
% GLABS = strcat( repmat("i",ng,1), ordr  , repmat("-",ng,1) , LOCI.GENE );


fields = fieldnames(ML);
MLname = fields{2};
MLtype = ML.(fields{2});

if strcmp(MLname,'ClassificationTree')


    % IMPORTANCE = predictorImportance(ML.ClassificationTree);
    IMPORTANCE = ML.ClassificationTree.predictorImportance;


    fh1 = figure('Units','normalized','Position',[.05 .04 .9 .9],...
        'Color','w','MenuBar','none');
    ax1 = axes('Position',[.06 .15 .9 .8],'Color','none');
    axes(ax1); ph1 = bar(rescale(IMPORTANCE,0,1));
    ax1.XTick = 1:ax1.XLim(2);
    ax1.XTickLabel = LOCI.GENE;
    ax1.XTickLabelRotation = 90;
    title('Predictor Importance');
    ylabel('Relative importance');


elseif strcmp(MLname,'ClassificationDiscriminant')

% fieldnames(MLtype)

    IMPORTANCE = MLtype.DeltaPredictor;

    WEIGHT = MLtype.W;


    fh1 = figure('Units','normalized','Position',[.05 .04 .9 .9],...
        'Color','w','MenuBar','none');
    ax1 = axes('Position',[.06 .15 .9 .8],'Color','none');
    axes(ax1); ph1 = bar(rescale(IMPORTANCE,0,1));
    ax1.XTick = 1:ax1.XLim(2);
    ax1.XTickLabel = LOCI.GENE;
    ax1.XTickLabelRotation = 90;
    title('Predictor Importance');
    ylabel('DeltaPredictor');




elseif strcmp(MLname,'GeneralizedLinearModel')

% fieldnames(MLtype)

    IMPORTANCE = MLtype.Coefficients.pValue(2:end);

    COEFF = MLtype.Coefficients.Estimate(2:end);


    fh1 = figure('Units','normalized','Position',[.05 .04 .9 .9],...
        'Color','w','MenuBar','none');
    ax1 = axes('Position',[.06 .15 .9 .8],'Color','none');
    axes(ax1); ph1 = bar(-log(IMPORTANCE));
    ax1.XTick = 1:ax1.XLim(2);
    ax1.XTickLabel = LOCI.GENE;
    ax1.XTickLabelRotation = 90;
    title('Predictor Importance');
    ylabel('-log(P)');


    fh2 = figure('Units','normalized','Position',[.05 .04 .9 .9],...
        'Color','w','MenuBar','none');
    ax2 = axes('Position',[.06 .15 .9 .8],'Color','none');
    axes(ax2); ph2 = bar(COEFF);
    ax2.XTick = 1:ax2.XLim(2);
    ax2.XTickLabel = LOCI.GENE;
    ax2.XTickLabelRotation = 90;
    title('Predictor Coefficient');
    ylabel('COEFF');


elseif strcmp(MLname,'ClassificationSVM')

% fieldnames(MLtype)

    IMPORTANCE = MLtype.Beta;

    if ~isempty(IMPORTANCE)

        fh1 = figure('Units','normalized','Position',[.05 .04 .9 .9],...
            'Color','w','MenuBar','none');
        ax1 = axes('Position',[.06 .15 .9 .8],'Color','none');
        axes(ax1); ph1 = bar(IMPORTANCE);
        ax1.XTick = 1:ax1.XLim(2);
        ax1.XTickLabel = LOCI.GENE;
        ax1.XTickLabelRotation = 90;
        title('Predictor Importance');
        ylabel('Beta');

    else

        IMPORTANCE = sum(MLtype.SupportVectors);

        fh1 = figure('Units','normalized','Position',[.05 .04 .9 .9],...
            'Color','w','MenuBar','none');
        ax1 = axes('Position',[.06 .15 .9 .8],'Color','none');
        axes(ax1); ph1 = bar(IMPORTANCE);
        ax1.XTick = 1:ax1.XLim(2);
        ax1.XTickLabel = LOCI.GENE;
        ax1.XTickLabelRotation = 90;
        title('Predictor Importance');
        ylabel('Support Vectors');
        

    end



elseif strcmp(MLname,'ClassificationKNN')


% fieldnames(MLtype)


%     IMPORTANCE = MLtype.NumNeighbors;
% 
%     fh1 = figure('Units','normalized','Position',[.05 .04 .9 .9],...
%         'Color','w','MenuBar','none');
%     ax1 = axes('Position',[.06 .15 .9 .8],'Color','none');
%     axes(ax1); ph1 = bar(IMPORTANCE);
%     ax1.XTick = 1:ax1.XLim(2);
%     ax1.XTickLabel = LOCI.GENE;
%     ax1.XTickLabelRotation = 90;
%     title('Predictor Importance');
%     ylabel('Support Vectors');


else


end


%%


% varargout = {FIT, PCTCORRECT, IMPORTANCE, LOCI.GENE, LOCI.CHRPOS};

end



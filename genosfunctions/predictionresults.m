function [varargout] = predictionresults(ML,TRX,TEX,LOCI)
%% predictionresults


TRAINX = TRX(:,10:end);
TRAINL = TRX(:,2);

TESTX = TEX(:,10:end);
TESTL = TEX(:,2);


FIT = ML.predictFcn(TESTX);

GUESSES = round(FIT) == TESTL;

PCTCORRECT = mean(GUESSES) .* 100;


[LABgi,~] = findgroups(TESTL);
[FITgi,~] = findgroups(round(FIT));



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



keyboard
%% ################   PREDICTOR IMPORTANCE   ################

% ng = numel(LOCI.GENE);
% ordr   = string((1:ng))';
% GHEAD = strcat( repmat("i",ng,1), ordr  , repmat("_",ng,1) , LOCI.GENE );
% GLABS = strcat( repmat("i",ng,1), ordr  , repmat("-",ng,1) , LOCI.GENE );

FIT    = ML.predictFcn(TRAINX);
fields = fieldnames(ML);
MLname = fields{2};
MLtype = ML.(fields{2});

Ntop = 50;


TRXCASEi = TRX(:,2)==1;
TRXCTRLi = TRX(:,2)~=1;
TRCASEnv = sum(TRX(TRXCASEi,10:end)>0);
TRCTRLnv = sum(TRX(TRXCTRLi,10:end)>0);
TROR = TRCASEnv ./ TRCTRLnv;
TRFISHOR = LOCI.TRFISHOR;


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
%%

        [KEY , KEYi] = sort(abs(IMPORTANCE),'descend');
        KEY = KEY(1:Ntop);
        KEYi = KEYi(1:Ntop);
        RISK = IMPORTANCE(KEYi)>0;

        clc; close all;
        fh1 = figure('Units','normalized','Position',[.05 .04 .9 .9],...
            'Color','w','MenuBar','none');
        ax1 = axes('Position',[.06 .15 .9 .8],'Color','none');
        axes(ax1); 
        ph1 = bar(KEY);
        ph1.FaceColor = 'flat';
        cm = lines(2);
        ph1.CData = cm(RISK+1,:);
        ax1.XTick = 1:Ntop;
        ax1.XTickLabel = LOCI.GENE(KEYi);
        ax1.XTickLabelRotation = 90;
        ax1.FontSize = 16;
        title('Predictor Importance');
        ylabel('Beta');
        ax1.YLim = [0 6.5];


        FISHOR = LOCI.TRFISHOR(KEYi);
        FISHOR = FISHOR(1:Ntop);

        FISHP  = LOCI.TRFISHP(KEYi);
        FISHP  = round(-log(FISHP(1:Ntop)),1);

        FCASEALT  = LOCI.TRCASEALT(KEYi);
        FCASEALT  = FCASEALT(1:Ntop);
        FCTRLALT  = LOCI.TRCTRLALT(KEYi);
        FCTRLALT  = FCTRLALT(1:Ntop);
        ALTS      = num2str([FCASEALT FCTRLALT]);
        text(1:Ntop,KEY,ALTS,'FontSize',12,'Rotation',90)

%%
    else


        SupportVectors      = MLtype.SupportVectors;
        SupportVectorLabels = MLtype.SupportVectorLabels;

        MX = SupportVectors .* SupportVectorLabels;

        IMPORTANCE = sum(MX);

        [KEY , KEYi] = sort(abs(IMPORTANCE),'descend');
        KEY = KEY(1:Ntop);
        KEYi = KEYi(1:Ntop);

        ODDS = TRFISHOR(KEYi);
        ODDS = ODDS(1:Ntop); ODDS(isnan(ODDS)|isinf(ODDS))=.5;
        %spf = sprintf('%0.1f\n',ODDS); ChrPosRefAlt = char(  string(spf)); 
        

        GENENAMES = [LOCI.GENE(KEYi) num2str(ODDS)];
        %strcat(char(GENENAMES),3)
        %strjoin(string(GENENAMES'))

        x(:,:,1) = cellstr(  GENENAMES(:,1)  );
        x(:,:,2) = cellstr(   num2str(ODDS)  );
        y = string(x);
        XLABS = join(y,3);


        clc; close all;
        fh1 = figure('Units','normalized','Position',[.05 .04 .9 .9],...
            'Color','w','MenuBar','none');
        ax1 = axes('Position',[.06 .15 .9 .8],'Color','none');
        axes(ax1); 
        ph1 = bar(KEY);
        ax1.XTick = 1:numel(KEY);
        %ax1.XTickLabel = LOCI.GENE(KEYi);    %'\fontsize{15} 
        ax1.XTickLabel = XLABS;    %'\fontsize{15} 
        ax1.XTickLabelRotation = 90;
        ax1.XTickLabel;
        ax1.FontSize = 14;
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



elseif strcmp(MLname,'ClassificationEnsemble')


    IMPORTANCE = MLtype;
    predictorImportance(MLtype)



else


end


%%


varargout = {FIT, PCTCORRECT, IMPORTANCE, LOCI.GENE, LOCI.CHRPOS};

end



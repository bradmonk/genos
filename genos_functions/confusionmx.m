function [CONMX, MU, XYAREA] = confusionmx(LABS,ACT,varargin)
%======================================================================
% BUILD CONFUSION MATRICES
%======================================================================
%  %CORRECT     NCASE_ACTUAL   NCTRL_ACTUAL   %CASE_ACTUAL
%  NCASE_PRED       TP             FP             PPV
%  NCTRL_PRED       FN             FN             NPV
%  %CASE_PRED      TPR            TNR              J
%------------------------------------------------------------------
% keyboard
%%

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

    % HIGH: PROPORTION POPULATION
    %HIPO_CASE = mean(ISHI_CASE);
    %HIPO_CTRL = mean(ISHI_CTRL);

    % HIGH: PROPORTION CORRECT
    HIMU_CASE = nanmean(ACT_CASE(ISHI_CASE) > 0);
    HIMU_CTRL = nanmean(ACT_CTRL(ISHI_CTRL) > 0);

    MU_ALL = nanmean([TF_CASE TF_CTRL]);
    HIMU_ALL = nanmean([(ACT_CASE(ISHI_CASE)>0) (ACT_CTRL(ISHI_CTRL)>0)]);


    MU = [MU_CASE MU_CTRL MU_ALL HIMU_CASE HIMU_CTRL HIMU_ALL];




    edges = (-.5:.02:.5);
    [CASE_NUM_LO,CASE_EDG_LO] = histcounts( (ACT_CASE(ACT_CASE>0)),edges);
    [CTRL_NUM_LO,CTRL_EDG_LO] = histcounts( (ACT_CTRL(ACT_CTRL>0).*-1),edges );
    [CASE_NUM_HI,CASE_EDG_HI] = histcounts( (ACT_CASE(ACT_CASE<0)),edges);
    [CTRL_NUM_HI,CTRL_EDG_HI] = histcounts( (ACT_CTRL(ACT_CTRL<0).*-1),edges );
    XAREA = mean([edges(1:end-1); edges(2:end)])';
    YAREA = [CASE_NUM_HI; CTRL_NUM_HI; CASE_NUM_LO; CTRL_NUM_LO]';

    XYAREA = [XAREA YAREA];


%%
% keyboard
%% GENERATE HISTOGRAM PLOTS

if nargin > 2
    close all
    fh02 = figure('Units','normalized','OuterPosition',[.05 .07 .60 .70],...
                  'Color','w','MenuBar','none');
    ax11 = axes('Position',[.06 .56 .4 .4],'Color','none');
    ax12 = axes('Position',[.56 .56 .4 .4],'Color','none');
    ax13 = axes('Position',[.06 .06 .4 .4],'Color','none');
    ax14 = axes('Position',[.56 .06 .4 .4],'Color','none');

    edges = (-.5:.02:.5);
    axes(ax11); histogram(ACT,(0:.02:1));      title('Neural Net Activation')
    axes(ax12); histogram(HACT,edges);         title('Classifier Zeroed Cutoff')
    axes(ax13); histogram(ACT_CASE,edges);     title('Case Activation')
    axes(ax14); histogram(ACT_CTRL.*-1,edges); title('Control Activation')
end


%% GENERATE AREA AND BAR PLOT
if nargin > 2
% close all

    fh1 = figure('Units','normalized','OuterPosition',[.01 .06 .55 .37],'Color','w');
    ax1 = axes('Position',[.07 .16 .61 .80],'Color','none');
    ax2 = axes('Position',[.78 .16 .21 .80],'Color','none');



    %---PLOT  AREAGRAM ---------------------
    axes(ax1) 
        area(XAREA,YAREA)
        legend({'Case Miss','Ctrl Miss',...
                'Case Hit','Ctrl Hit'},...
                'Location','best');
        ylabel('Count')
        ax1.FontSize=16;
        %title('Classifier Performance Area');
    %----------------------------------------------------------------------


    %---BAR GRAPH ---------------------
    axes(ax2)
    bar(([ HIMU_CASE , HIMU_CTRL , HIMU_ALL ].*100),.5, 'FaceColor',[.31 .31 .31]); 
        hold on; 
    bar(([ MU_CASE, MU_CTRL, MU_ALL ].*100),.20,'FaceColor',[.95 .85 .50]);
    
    %---BAR GRAPH FORMATTING ---------------------
    grid on; ylabel('Pct. Correct')
    legend({'Top 25%','All'},'Location','Northwest','NumColumns',2)
    ax2.YLim = [0 116]; 
    ax2.YTick = [0 25 50 75 100]; 
    ax2.XTickLabels = {'CASE','CTRL','ALL'};
    ax2.XTickLabelRotation = 33;
    ax2.FontSize=16;
    %title('Performance Summary')
    %----------------------------------------------------------------------
    pause(.2)



    %ax5=subplot(3,7,20:21); 
    %    bar([h3.Values; h4.Values; h1.Values; h2.Values]','stacked')
    %XYAREA = [XAREA YAREA];
    %varargout = {XYAREA};
end
%%


end



% sum(HO_DL(:,1))
% sum(~HO_DL(:,1))
% sum(round(HO_ACT(1,:)))
% sum(round(HO_ACT(2,:)))
% clc;
% disp(CONMX_TR)
% disp(' ')
% disp(CONMX_HO)
% tabulate(categorical(HO_DL(:,1)))
% plotconfusion(DL',ACT)
% output = plotconfusion(HO_DL',HO_ACT);
% HODL  = HO_DL'; HODL = HODL(1,:); HODL(HODL==0) = 2;
% HOACT = round(HO_ACT(1,:)); HOACT(HOACT==0) = 2;
% C = confusionmat(HODL,HOACT)';
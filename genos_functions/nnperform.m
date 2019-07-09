function [CONMX, MU, MU_CUT, XYAREA] = nnperform(YL,MX,PVX,net,varargin)
%======================================================================
% Parse Function Inputs
%======================================================================

%net.userdata.note

p = inputParser;
defaultCut = .5;
defaultPlot = 0;
addOptional(p,'cut',defaultCut,@isnumeric)
addOptional(p,'plt',defaultPlot,@isnumeric)
parse(p,varargin{:})


NETDAT.CUT = p.Results.cut;
NETDAT.PLT = p.Results.plt;


NETDAT.PVX = PVX;


%======================================================================
%% GET NEURAL NET OUTPUTS FOR TRAINING AND HOLDOUT SET
%======================================================================


ACT = net(MX);



%==========================================================================
%% DETERMINE CLASSIFIER CUTOFF BASED ON POPULATION PERCENTAGES
%==========================================================================
clearvars -except YL MX net VIN ACT NETDAT


CACOn = numel(YL);
CTRLn = sum(~YL);
CASEn = sum(YL);


% ACTUAL PROPORTION OF CASES & CTRLS
actualCASEp = mean(YL);
actualCTRLp = mean(~YL);


% GUESSED PROPORTION OF CASES & CTRLS AT CUT=0.5
guessCASEp = mean(ACT >  .5);
guessCTRLp = mean(ACT <= .5);


% GET CUTOFF TO MAKE ACTUAL & GUESSED PROPORTIONS MORE SIMILAR
CASEcut = quantile( ACT, actualCASEp );
CTRLcut = quantile( ACT, actualCTRLp );
cut = mean([CASEcut CTRLcut]);


% GUESSED PROPORTION OF CASES & CTRLS AT CUT=0.5
guessCASEScut = mean(ACT >  cut);
guessCTRLScut = mean(ACT <= cut);




clc; disp(' ');
disp(net.userdata.note)
fprintf('\n CACO N: %.0f \n CTRL N: %.0f \n CASE N: %.0f \n', CACOn,CTRLn,CASEn)

fprintf('\n WHEN CUTOFF IS %.2f ...\n',0.5);
fprintf(' CTRL P(ACTUAL,GUESS): %.2f, %.2f \n', actualCTRLp, guessCTRLp)
fprintf(' CASE P(ACTUAL,GUESS): %.2f, %.2f \n', actualCASEp, guessCASEp)


fprintf('\n WHEN CUTOFF IS %.2f ...\n',cut);
fprintf(' CTRL P(ACTUAL,GUESS): %.2f, %.2f \n', actualCTRLp, guessCTRLScut)
fprintf(' CASE P(ACTUAL,GUESS): %.2f, %.2f \n', actualCASEp, guessCASEScut)



NETDAT.Qcut = cut;


%==========================================================================
%% MAKE CONFUSION MATRIX
%==========================================================================
%  %CORRECT     NCASE_ACTUAL   NCTRL_ACTUAL   %CASE_ACTUAL
%  NCASE_PRED       TP             FP             PPV
%  NCTRL_PRED       FN             FN             NPV
%  %CASE_PRED      TPR            TNR              J
%--------------------------------------------------------------------------
clearvars -except YL MX net VIN ACT cut NETDAT


if NETDAT.CUT == 1

    OUT = ACT + (.5-NETDAT.Qcut);
    OUT(OUT>1)=1;
    OUT(OUT<0)=0;

else
    OUT = ACT;
end


[ERR, CMX,~, CMR] = confusion(YL,OUT);
COR = 1-ERR;
CMX = CMX';


CONMX = nan(4,4);

CONMX(2:3,2:3) = CMX;

CONMX(2:3,1) = [CMX(1,1)+CMX(1,2);
                CMX(2,1)+CMX(2,2)];

CONMX(1,2:3) = [CMX(1,1)+CMX(2,1),...
                CMX(1,2)+CMX(2,2)];


TPR = CONMX(2,2) / CONMX(1,2);
TNR = CONMX(3,3) / CONMX(1,3);
PPV = CONMX(2,2) / CONMX(2,1);
NPV = CONMX(3,3) / CONMX(3,1);

PACASE = CONMX(1,2) / (CONMX(1,2) + CONMX(1,3));
PPCASE = CONMX(2,1) / (CONMX(2,1) + CONMX(3,1));

CONMX(1,4) = PACASE;
CONMX(4,1) = PPCASE;
CONMX(4,2) = TPR;
CONMX(4,3) = TNR;
CONMX(2,4) = PPV;
CONMX(3,4) = NPV;
CONMX(1,1) = COR;
CONMX(4,4) = TPR + TNR - 1;

NETDAT.CONMX = CONMX;

disp(CONMX);



% CHECK IF PREDICTION COUNTS MATCH
%all(sum(round(ACTt),2) == CONMX(2:3,1))

% CHECK IF ACTUAL LABEL COUNTS MATCH
%all(sum(LABSt,2)' == CONMX(1,2:3))
%==========================================================================
%% DETERMINE NEURAL NET PERFORMANCE FOR TRAINING DATA
%==========================================================================
clearvars -except YL MX net VIN ACT cut NETDAT


if NETDAT.CUT == 1
    CUTOFF = NETDAT.Qcut;
else
    CUTOFF = .5;
end


OUT = ACT;

% SPLIT PREDICTIONS FOR CASES AND CONTROLS
ACT_CASE = OUT(YL==1);
ACT_CTRL = OUT(YL~=1);


% PROPORTION CORRECT AT CUTOFF
TF_CASE = ACT_CASE >  CUTOFF; 
TF_CTRL = ACT_CTRL <= CUTOFF; 

MU_CASE = nanmean(TF_CASE);
MU_CTRL = nanmean(TF_CTRL);
MU_CACO = nanmean([TF_CASE TF_CTRL]);


T_CASE_ACT = ACT_CASE( TF_CASE );
F_CASE_ACT = ACT_CASE(~TF_CASE );
T_CTRL_ACT = ACT_CTRL( TF_CTRL );
F_CTRL_ACT = ACT_CTRL(~TF_CTRL );



% PROPORTION CORRECT AT QCUT
TF_CASE_cut = ACT_CASE >  NETDAT.Qcut;
TF_CTRL_cut = ACT_CTRL <= NETDAT.Qcut;

MU_CASE_cut = nanmean(TF_CASE_cut);
MU_CTRL_cut = nanmean(TF_CTRL_cut);
MU_CACO_cut = nanmean([TF_CASE_cut TF_CTRL_cut]);


T_CASE_ACT_cut = ACT_CASE( TF_CASE_cut );
F_CASE_ACT_cut = ACT_CASE(~TF_CASE_cut );
T_CTRL_ACT_cut = ACT_CTRL( TF_CTRL_cut );
F_CTRL_ACT_cut = ACT_CTRL(~TF_CTRL_cut );



fprintf('\n ABSOLUTE CUTOFF [%.2f] (CACO|CTRL|CASE): %.2f | %.2f | %.2f   ',...
    CUTOFF, MU_CACO, MU_CTRL, MU_CASE)
fprintf('\n POP-FREQ CUTOFF [%.2f] (CACO|CTRL|CASE): %.2f | %.2f | %.2f \n',...
    NETDAT.Qcut, MU_CACO_cut, MU_CTRL_cut, MU_CASE_cut)






% USE QUANTILE FUNCTION TO ESTABLISH A HIGH-CONFIDENCE CUTOFF
Q = quantile(OUT,[.25,.75]);
LOQ = Q(1); HIQ = Q(2);

% GET PROPORTION CORRECT FOR HIGH CONFIDENCE OUTPUTS
TF_CASE_HI = OUT( (YL==1) & ((OUT < LOQ)|(OUT > HIQ))) >   CUTOFF;
TF_CTRL_HI = OUT( (YL~=1) & ((OUT < LOQ)|(OUT > HIQ))) <=  CUTOFF;

MU_CASE_HI = nanmean(TF_CASE_HI);
MU_CTRL_HI = nanmean(TF_CTRL_HI);
MU_CACO_HI = nanmean([TF_CTRL_HI TF_CASE_HI]);


% GET PROPORTION CORRECT FOR HIGH CONFIDENCE OUTPUTS
TF_CASE_HI = OUT( (YL==1) & ((OUT < LOQ)|(OUT > HIQ))) >   NETDAT.Qcut;
TF_CTRL_HI = OUT( (YL~=1) & ((OUT < LOQ)|(OUT > HIQ))) <=  NETDAT.Qcut;

MU_CASE_HI_cut = nanmean(TF_CASE_HI);
MU_CTRL_HI_cut = nanmean(TF_CTRL_HI);
MU_CACO_HI_cut = nanmean([TF_CTRL_HI TF_CASE_HI]);


fprintf('\n HICONF [%.2f | %.2f] ABSOLUTE CUTOFF [%.2f] (CACO|CTRL|CASE): %.2f | %.2f | %.2f   ',...
    LOQ, HIQ, CUTOFF, MU_CACO_HI, MU_CTRL_HI, MU_CASE_HI)
fprintf('\n HICONF [%.2f | %.2f] POP-FREQ CUTOFF [%.2f] (CACO|CTRL|CASE): %.2f | %.2f | %.2f \n',...
    LOQ, HIQ, NETDAT.Qcut, MU_CACO_HI_cut, MU_CTRL_HI_cut, MU_CASE_HI_cut)


MU     = [MU_CASE MU_CTRL MU_CACO MU_CASE_HI MU_CTRL_HI MU_CACO_HI];
MU_CUT = [MU_CASE_cut MU_CTRL_cut MU_CACO_cut MU_CASE_HI_cut MU_CTRL_HI_cut MU_CACO_HI_cut];





% GET NEURAL NET OUTPUT AREA HISTOGRAM DATA
edges = (-.5:.02:.5);
CASE_NUM_LO = histcounts( F_CASE_ACT - .5 , edges );
CTRL_NUM_LO = histcounts( F_CTRL_ACT - .5 , edges );
CASE_NUM_HI = histcounts( T_CASE_ACT - .5 , edges );
CTRL_NUM_HI = histcounts( T_CTRL_ACT - .5 , edges );
XAREA = mean([edges(1:end-1); edges(2:end)])';
YAREA = [CASE_NUM_LO; CTRL_NUM_LO; CASE_NUM_HI; CTRL_NUM_HI]';
XYAREA = [XAREA YAREA];




NETDAT.MU     = MU;
NETDAT.MU_CUT = MU_CUT;
NETDAT.XAREA  = XAREA;
NETDAT.YAREA  = YAREA;
NETDAT.XYAREA = XYAREA;




%==========================================================================
%% PLOT DATA
%==========================================================================
clearvars -except YL MX net VIN ACT cut NETDAT
if NETDAT.PLT==1

    nnperfplots(YL,MX,ACT,NETDAT)

end





%==========================================================================
%% EXPORT VARIABLES
%==========================================================================

clearvars -except YL MX net VIN ACT cut NETDAT


CONMX  = NETDAT.CONMX;
MU     = NETDAT.MU;
MU_CUT = NETDAT.MU_CUT;
XYAREA = NETDAT.XYAREA;





% close all;
% fh1 = figure('Units','normalized','Position',[.02 .04 .95 .85],'Color','w');
% ax1 = axes('Position',[.05 .53 .93 .45],'Color','none');
% ax2 = axes('Position',[.05 .04 .93 .45],'Color','none');
% axes(ax1); bar( ([sum(VXt==-1);sum(VXt==0);sum(VXt==2);sum(VXt==3)]') ,'stacked');axis tight
% axes(ax2); bar( ([sum(DXt==-1);sum(DXt==0);sum(DXt==2);sum(DXt==3)]') ,'stacked');axis tight

end
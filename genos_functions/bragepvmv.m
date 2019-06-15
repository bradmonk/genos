function [PVMX,varargout] = bragepvmv(PV,varargin)
%%  AGEZ & BRAGEZ EXPLAINED
% 
%{

When training a classifier to recognize AD, it should
*weight heavily*...

- YOUNG CASES 
- OLDER CONTROLS

- HIGH BRAAK CASES
- LOW BRAAK CONTROLS


We can design a variable called 'BRAGE' that is based on both 
Age and a BRAD function...

    BRAGE = AGE - BRAD

...where 

    BRAD = tansig(-1.5:.5:2).*10

...the worse their 'BRAAK' value the more their age is reduced:

    BRAGE(BRAAK==0)   : AGE + 9.1
    BRAGE(BRAAK==1)   : AGE + 7.6
    BRAGE(BRAAK==2)   : AGE + 4.6
    BRAGE(BRAAK==NaN) : AGE - 0.0
    BRAGE(BRAAK==3)   : AGE - 4.6
    BRAGE(BRAAK==4)   : AGE - 7.6
    BRAGE(BRAAK==5)   : AGE - 9.1
    BRAGE(BRAAK==6)   : AGE - 9.6






So for example if someone is 80 years old and has a BRAAK==0 their
BRAGE equals 89.1; if they had a BRAAK==6 their BRAGE equals 70.4

From there a Z-score is computed for both AGE & BRAGE.

!!IMPORTANTLY!! YOU WILL PROBABLY WANT TO FLIP THE *CASE* AGEz & BRAGEz
BEFORE USING (the sign for CASES is flipped to make 
negative-AGEz and negative-BRAGEz CASES high positive values 
(highly weighted).


%}
% 
%--------------------------------------------------------------------------
% COLS:  PVMX( 1 , 2, 3 , 4 ,  5 , 6 ,  7  ,  8  ,  9    ...)
% START: PVMX(SRR,AD,COH,AGE,APOE,SEX,BRAAK,BRAAK,BRAAK.....)
% END:   PVMX(SRR,AD,COH,AGE,APOE,BRAD,BRAAK,AGEz,BRAGEz....)
% 
%%

PVMX = PV;

if nargin == 2
    MakeCaseBrageNegative = 1;
else
    MakeCaseBrageNegative = 0;
end




%% MAKE U55 CASES EQUAL TO 55; MAKE U60 CONTROLS EQUAL TO 60
% COLS:  PVMX( 1 , 2, 3 , 4 ,  5 , 6 ,  7  ,  8  ,  9    ...)
% START: PVMX(SRR,AD,COH,AGE,APOE,SEX,BRAAK,BRAAK,BRAAK.....)
% END:   PVMX(SRR,AD,COH,AGE,APOE,BRAD,BRAAK,AGEz,BRAGEz....)

CA = (PVMX(:,2) == 1);
CO = (PVMX(:,2) ~= 1);

U55 = (PVMX(:,4) < 55);
U60 = (PVMX(:,4) < 60);


PVMX(CA&U55,4) = 55;
PVMX(CO&U60,4) = 60;



%% COMPUTE THE Z-SCORE OF AGE AND OVERWRITE COL *8*
% COLS:  PVMX( 1 , 2, 3 , 4 ,  5 , 6 ,  7  ,  8  ,  9    ...)
% START: PVMX(SRR,AD,COH,AGE,APOE,SEX,BRAAK,BRAAK,BRAAK.....)
% END:   PVMX(SRR,AD,COH,AGE,APOE,BRAD,BRAAK,AGEz,BRAGEz....)


PVMX(CA,8) = zscore(PVMX(CA,4));
PVMX(CO,8) = zscore(PVMX(CO,4));




%==========================================================================
%% CREATE "BRAGE" (BRAAK-BASED-AGE):  BRAGE = BRAAK-AGE
%==========================================================================
% COLS:  PVMX( 1 , 2, 3 , 4 ,  5 , 6 ,  7  ,  8  ,  9    ...)
% START: PVMX(SRR,AD,COH,AGE,APOE,SEX,BRAAK,BRAAK,BRAAK.....)
% END:   PVMX(SRR,AD,COH,AGE,APOE,BRAD,BRAAK,AGEz,BRAGEz....)
% [BRAX] = bragepvmv(PVTR);



% SET COL*9* BRAAK NaN TO -10
%------------------------------------------------
ISNA = isnan(PVMX(:,9));
PVMX(CA&ISNA,9) = -10;
PVMX(CO&ISNA,9) = -10;


BRADVALS = round(tansig(-1.5:.5:2).*10,1);

% BRAD = BRAAK-AGE SCORING SYSTEM
%------------------------------------------------
CAi = CA&(PVMX(:,9)==[0:2 -10 4:7]);
COi = CO&(PVMX(:,9)==[0:2 -10 4:7]);

PVMX( CAi(:,1) ,9)  = BRADVALS(1);
PVMX( COi(:,1) ,9)  = BRADVALS(1);

PVMX( CAi(:,2) ,9)  = BRADVALS(2);
PVMX( COi(:,2) ,9)  = BRADVALS(2);

PVMX( CAi(:,3) ,9)  = BRADVALS(3);
PVMX( COi(:,3) ,9)  = BRADVALS(3);

PVMX( CAi(:,4) ,9)  = BRADVALS(4);
PVMX( COi(:,4) ,9)  = BRADVALS(4);

PVMX( CAi(:,5) ,9)  = BRADVALS(5);
PVMX( COi(:,5) ,9)  = BRADVALS(5);

PVMX( CAi(:,6) ,9)  = BRADVALS(6);
PVMX( COi(:,6) ,9)  = BRADVALS(6);

PVMX( CAi(:,7) ,9)  = BRADVALS(7);
PVMX( COi(:,7) ,9)  = BRADVALS(7);

PVMX( CAi(:,8) ,9)  = BRADVALS(8);
PVMX( COi(:,8) ,9)  = BRADVALS(8);









% BRAGE = AGE - BRAD
%------------------------------------------------
% COLS:  PVMX( 1 , 2, 3 , 4 ,  5 , 6 ,  7  ,  8  ,  9    ...)
% START: PVMX(SRR,AD,COH,AGE,APOE,SEX,BRAAK,BRAAK,BRAAK.....)
% END:   PVMX(SRR,AD,COH,AGE,APOE,BRAD,BRAAK,AGEz,BRAGEz....)
% [BRAX] = bragepvmv(PVTR);


% MAKE COL6 BRAD ARRAY
PVMX(:,6) = PVMX(:,9);


% MAKE COL9 BRAGE
PVMX(:,9) = PVMX(:,4) - PVMX(:,6);


% MAKE COL6 BRAGE
PVMX(:,6) = PVMX(:,9);


% MAKE COL9 BRAGEz
PVMX(CA,9) = zscore(PVMX(CA,9));
PVMX(CO,9) = zscore(PVMX(CO,9));

% MAKE COL8 AGEz
PVMX(CA,8) = zscore(PVMX(CA,8));
PVMX(CO,8) = zscore(PVMX(CO,8));




% APPLY BRAGE WEIGHT ALTERATION TO SNPS
%------------------------------------------------
PVAX = PVMX;
CA = (PVAX(:,2) == 1);
PVAX(CA,9) = PVAX(CA,9).*-1;
SNPS = PVAX(:,10:end) > 0;
PVAX(:,10:end) = PVAX(:,10:end) + (PVAX(:,9).*SNPS);
varargout = {PVAX};
%------------------------------------------------



%------------------------------------------------
if MakeCaseBrageNegative
PVMX(:,9) = BRAGE;
end
%------------------------------------------------




%clearvars -except ADSP LOCI CASE CTRL PHEN USNP PHECASE PHECTRL
%% PLOT HISTOGRAM OF AGE AND AGE-TRANSFORMATIONS
%------------------------------------------------
% END:   PVMX(SRR,AD,COH,AGE,APOE,BRAGE,BRAAK,AGEz,BRAGEz....)
% BRAX = bragepvmv(PVTR);
%{
close all;
fh01 = figure('Units','normalized','OuterPosition',[.01 .05 .95 .90],...
              'Color','w','MenuBar','none');
ax01 = axes('Position',[.06 .56 .4 .4],'Color','none');
ax02 = axes('Position',[.56 .56 .4 .4],'Color','none');
ax03 = axes('Position',[.06 .06 .4 .4],'Color','none');
ax04 = axes('Position',[.56 .06 .4 .4],'Color','none');


rbins=50.5:1:90.5;
zbins=-4:.25:4;

axes(ax01); histogram(PVMX(CA,4),rbins);
hold on;    histogram(PVMX(CO,4),rbins); xlabel('CASE vs. CTRL AGE')
axes(ax02); histogram(PVMX(CA,8),zbins);
hold on;    histogram(PVMX(CO,8),zbins); xlabel('CASE vs. CTRL AGE Z-SCORE')
axes(ax03); histogram(PVMX(CA,6),rbins);
hold on;    histogram(PVMX(CO,6),rbins); xlabel('CASE vs. CTRL BRAGE')
axes(ax04); histogram(PVMX(CA,9),zbins);
hold on;    histogram(PVMX(CO,9),zbins); xlabel('CASE vs. CTRL BRAGE Z-SCORE')

%}


%clearvars -except ADSP LOCI CASE CTRL PHEN USNP PHECASE PHECTRL
%%
end





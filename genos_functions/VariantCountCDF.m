%###############################################################
%%  USE PER-PERSON VARIANT COUNTS TO IDENTIFY OUTLIERS
%###############################################################

PHECASE = PHEN(PHEN.AD==1,:);
PHECTRL = PHEN(PHEN.AD~=1,:);

PHECASE(PHECASE.TOTvars < 10000, : ) = [];
PHECTRL(PHECTRL.TOTvars < 10000, : ) = [];


Qlo = quantile(PHECASE.TOTvars,.01);
Qhi = quantile(PHECASE.TOTvars,.99);
PHECA  = PHECASE(((PHECASE.TOTvars > Qlo) & (PHECASE.TOTvars < Qhi)),:);


Qlo = quantile(PHETRCTRL.TOTvars,.01);
Qhi = quantile(PHETRCTRL.TOTvars,.99);
PHECO  = PHETRCTRL(((PHETRCTRL.TOTvars > Qlo) & (PHETRCTRL.TOTvars < Qhi)),:);



% clc; close all
% fh1=figure;
% ax1=subplot(2,1,1); ph1=histogram(PHECA.TOTvars,20,'Normalization','cdf'); title('CASE')
% ax2=subplot(2,1,2); ph2=histogram(PHECO.TOTvars,20,'Normalization','cdf'); title('CTRL')



clc; close all
fh1=figure;
ax1=subplot(2,1,1); ph1=histogram(PHECASE.TOTvars,20,'Normalization','cdf'); title('CASE')
ax2=subplot(2,1,2); ph2=histogram(PHECTRL.TOTvars,20,'Normalization','cdf'); title('CTRL')


bins = 17000:100:20000;
ph1.BinEdges = bins;
ph2.BinEdges = bins;
[bins(1:end-1)' ph1.Values' ph2.Values']






% close all;
% cohcounts(TRCASE,TRCTRL,TECASE,TECTRL,2)

% clearvars -except ADSP LOCI CASE CTRL PHEN USNP...
% TRCASE TRCTRL TECASE TECTRL

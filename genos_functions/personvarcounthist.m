function [] = personvarcounthist(PPVCOUNTS, COHNUM, varargin)
%% COID CONSORT	STUDY_NAME  STUDY   FAM PEOPLE  ENR     CASE	CTRL	TOT    %CA %CO  EQL	CULL	NOTES
% 1     DGC     Adult Chng  ACT     0	1277	0       323     945     1268	25	75	323	-622	KEEP
% 2     ADGC	AD Centers  ADC     0	3263	0       2438	817     3255	75	25	817	1621	KEEP
% 3     CHARGE	Athrosclro  ARIC	0	57      0       39      18      57      68	32	18	21      REMOVE
% 4     CHARGE	Aus Stroke	ASPS	0	126     0       121     5       126     96	4	5	116     REMOVE
% 5     ADGC	Chic Aging  CHAP	0	231     0       27      204     231     12	88	27	-177	REMOVE
% 6     CHARGE	CardioHlth  CHS     0	840     0       250     583     833     30	70	250	-333	KEEP
% 7     ADGC	Hispanic S  CUHS	361	0       331     160     171     331     48	52	160	-11     HISPANIC
% 8     CHARGE	Erasmus Er  ERF     30	45      0       45      0       45      100	0	0	45      REMOVE
% 9     CHARGE	Framingham  FHS     0	581     0       157     424     581     27	73	157	-267	KEEP
% 10	ADGC	Gene Diffs  GDF     0	207     0       111     96      207     54	46	96	15      KEEP
% 11	ADGC	NIA - LOAD  LOAD	80  113     363     367     109     476     77	23	109	258     KEEP
% 12	ADGC	Aging Proj  MAP     0	415     0       138     277     415     33	67	138	-139	KEEP
% 13	ADGC	Mayo Clini  MAYO	0	349     0       250     99      349     72	28	99	151     KEEP
% 14	ADGC	Miami Univ  MIA     61	181     19      186     14      200     93	7	14	172     REMOVE
% 15	ADGC	AD Genetic  MIR     0	284     47      316     15      331     95	5	15	301     REMOVE
% 16	ADGC	Mayo cl PD  MPD     0	20      0       0       20      20      0	100	0	-20     REMOVE
% 17	ADGC	NationC AD  NCRD	18	108     52      160     0       160     100	0	0	160     REMOVE
% 18	ADGC	Wash Unive  RAS     0	0       46      46      0       46      100	0	0	46      REMOVE
% 19	ADGC	Relig Ordr  ROS     0	351     0       154     197     351     44	56	154	-43     KEEP
% 20	CHARGE	RotterdamS  RS      0	1089	0       276     813     1089	25	75	276	-537	KEEP
% 21	ADGC	Texas AD S  TARC	0	144     0       132     12      144     92	8	12	120     REMOVE
% 22	ADGC	Un Toronto  TOR     0	0       9       9       0       9       100	0	0	9       REMOVE
% 23	ADGC	Vanderbilt  VAN     6	235     1       210     26      236     89	11	26	184     REMOVE
% 24	ADGC	WashNY Age  WCAP	0	147     3       34      116     150     23	77	34	-82     REMOVE
% 25	ADGC 	Univ Penns  UPN     40 	0       3       0       0       0       0	0	0	0       REMOVE


if nargin > 2
    lo = varargin{1}(1);
    hi = varargin{1}(2);
else
    lo = 0;
    hi = inf;
end

%%


PPVCOUNTS = PPVCOUNTS( (PPVCOUNTS>lo) & (PPVCOUNTS<hi) , :);


clc;
studs = {
'ACT'
'ADC'
'ARIC'
'ASPS'
'CHAP'
'CHS'
'CUHS'
'ERF'
'FHS'
'GDF'
'LOAD'
'MAP'
'MAYO'
'MIA'
'MIR'
'MPD'
'NCRD'
'RAS'
'ROS'
'RS'
'TARC'
'TOR'
'VAN'
'WCAP'
};

PPC.i = {};
PPC.v = {};
PPC.s = [];

for nn = 1:24

    PPC.i{nn} = COHNUM==nn;
    PPC.v{nn} = PPVCOUNTS( PPC.i{nn} );
    PPC.s(nn) = sum( PPC.v{nn} );

end




close all
fh0 = figure('Units','normalized','OuterPosition',[.1 .05 .7 .92],'Color','w');
ax0 = axes('Position',[.06 .06 .9 .9],'Color','none');
leg=[];
for nn = 1:24
    if any(PPC.v{nn})
    histogram(PPC.v{nn}, 40,'DisplayStyle','stairs',...
    'LineWidth',3,'Normalization','count');
    hold on; leg(end+1)=nn;
    end
end
legend(studs{leg},'Location','northwest')
ylabel('count')


fh1 = figure('Units','normalized','OuterPosition',[.1 .05 .7 .92],'Color','w');
ax1 = axes('Position',[.06 .06 .9 .9],'Color','none');
leg=[];
for nn = 1:24
    if any(PPC.v{nn})
    histogram(PPC.v{nn}, 40,'DisplayStyle','stairs',...
    'LineWidth',3,'Normalization','probability');
    hold on; leg(end+1)=nn;
    end
end
legend(studs{leg},'Location','northwest')
ylabel('probability')





fh1 = figure('Units','normalized','OuterPosition',[.1 .05 .7 .92],'Color','w');
ax1 = axes('Position',[.06 .06 .9 .9],'Color','none');
leg=[];
for nn = 1:24
    if any(PPC.v{nn})
    histogram(PPC.v{nn}, 40,'DisplayStyle','stairs',...
    'LineWidth',3,'Normalization','countdensity');
    hold on; leg(end+1)=nn;
    end
end
legend(studs{leg},'Location','northwest')
ylabel('countdensity')



fh1 = figure('Units','normalized','OuterPosition',[.1 .05 .7 .92],'Color','w');
ax1 = axes('Position',[.06 .06 .9 .9],'Color','none');
leg=[];
for nn = 1:24
    if any(PPC.v{nn})
    histogram(PPC.v{nn}, 40,'DisplayStyle','stairs',...
    'LineWidth',3,'Normalization','pdf');
    hold on; leg(end+1)=nn;
    end
end
legend(studs{leg},'Location','northwest')
ylabel('pdf')




fh1 = figure('Units','normalized','OuterPosition',[.1 .05 .7 .92],'Color','w');
ax1 = axes('Position',[.06 .06 .9 .9],'Color','none');
leg=[];
for nn = 1:24
    if any(PPC.v{nn})
    histogram(PPC.v{nn}, 40,'DisplayStyle','stairs',...
    'LineWidth',3,'Normalization','cdf');
    hold on; leg(end+1)=nn;
    end
end
legend(studs{leg},'Location','northwest')
ylabel('cdf')










end

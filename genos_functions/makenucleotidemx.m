function [NUCNUM, NUCTXT] = makenucleotidemx(ADMX,CASE,CTRL,USNP,PHE,varargin)
%{

SAMPLE.ped    (NOTE: REMOVE HEADER ROW IN ACTUAL .ped FILE)
--------------------------
PED       ID      DAD  MOM  SEX   AD   VAR1   VAR1   VAR3   VAR4...  
200001    500001    0    0    1    1    1 2    1 2    1 2    1 2...
200001    500002    0    0    1    1    3 3    3 2    1 2    1 2...
200001    500003    0    0    1    1    4 2    1 2    4 2    1 2...
200001    500001    0    0    1    1    1 2    1 2    1 2    1 2...
200001    500002    0    0    1    1    3 3    3 2    1 2    1 2...
200001    500003    0    0    1    1    4 2    1 2    4 2    1 2...
200001    500001    0    0    1    1    1 2    1 2    1 2    1 2...
200001    500002    0    0    1    1    3 3    3 2    1 2    1 2...
200001    500003    0    0    1    1    4 2    1 2    4 2    1 2...




SAMPLE.info    (NOTE: REMOVE HEADER ROW IN ACTUAL .ped FILE)
--------------------------
CHR      POS
chr01    274044
chr01    274149
chr01    274333
chr01    274718
chr01    275141
chr02    274333
chr02    274718
chr02    275141
chr02    279049
chr02    281011


variants with p-values < 1e-5


SAMPLE.ped  (Linkage Format) DEFINITIONS
----------------------------------------------------

PED: Pedigree Name A unique alphanumeric identifier for this individual's
family. Unrelated individuals should not share a pedigree name.

ID: Individual ID An alphanumeric identifier for this individual. Should be
unique within his family (see above).

DAD: Father's ID Identifier corresponding to father's individual ID or "0"
if unknown father. Note that if a father ID is specified, the father must
also appear in the file.

MOM: Mother's ID Identifier corresponding to mother's individual ID or "0"
if unknown mother Note that if a mother ID is specified, the mother must
also appear in the file.

SEX: Sex Individual's gender (1=MALE, 2=FEMALE).

AD: Case/Control Affection status Affection status to be used for
association tests (0=UNKNOWN, 1=UNAFFECTED,2=AFFECTED).

GENO: Marker genotypes Each marker is represented by two columns (one for
each allele, separated by a space) and coded either ACGT or 1-4 where: 1=A,
2=C, 3=G, T=4, Unknown=0


It is also worth noting that this format can be used with non-family based
data. Simply use a dummy value for the pedigree name (1, 2, 3...) and fill
in zeroes for father and mother ID. It is important that the "dummy" value
for the ped name be unique for each individual. Affection status can be
used to designate cases vs. controls (2 and 1, respectively).

Files should also follow the following guidelines:

Families should be listed consecutively within the file (i.e. all the lines
with the same pedigree ID should be adjacent) If an individual has a
nonzero parent, the parent should be included in the file on his own line.






SAMPLE.info  (Marker Information File) DEFINITIONS
----------------------------------------------------

The marker info file is two columns, marker name and position. The
positions can be either absolute chromosomal coordinates or relative
positions. It might look something like this:

marker01 190299 marker02 190950 marker03 191287 An optional third column
can be included in the info file to make additional notes for specific
SNPs. SNPs with additional information are highlighted in green on the LD
display. For instance, you could make note that the first SNP is a coding
variant as follows:

marker01 190299 CODING_SNP marker02 190950 marker03 191287




%}


% keyboard

if nargin > 5
    REFUNKALT = varargin{1};
else
              %refref unkunk refalt altalt
    REFUNKALT = [-1     -1     1     3];
end


CACOID = PHE.SRR;

vMX  = zeros( size(CACOID,1) , size(ADMX,1) );


for nn = 1:size(CASE,1)

    CASES = CASE{nn};
    CTRLS = CTRL{nn};

    CACO = [CASES; CTRLS];
    
    UNCAS = USNP{nn};




    if any(UNCAS)
        [~,Uj] = ismember(CACOID , UNCAS );
        vMX(Uj>0,nn) = 1;
    end
    if any(CACO)
        CACOSRR = CACO(:,1);
        CACOHH  = CACO(:,2);
    
        [~,Aj] = ismember(CACOID , CACOSRR );
        [~,Uj] = ismember(CACOID , UNCAS );

        Af = Aj(Aj>0);

        vMX(Aj>0,nn) = (CACOHH(Af)+1); %HETALT:2  HOMALT:3
        vMX(Uj>0,nn) = 1;              %UNCALL:1  HOMREF:0
    end


end

cacoMX = vMX;

cacoMX(vMX==0) =  REFUNKALT(1);     % HOMREF: -1
cacoMX(vMX==1) =  REFUNKALT(2);     % UNCALL: -1 or 0
cacoMX(vMX==2) =  REFUNKALT(3);     % HETALT:  1
cacoMX(vMX==3) =  REFUNKALT(4);     % HOMALT:  3



%------------------------------------------------------------
%----- ANNOTATE THE MATRIX WITH USER INFO
% ADD 6 COLUMNS AND 1 ROW TO THE MATRIX
ADNN = padarray(cacoMX,[1 6],0,'pre');

sz = size(ADNN(2:end,1),1);

% ADNN(2:end , 1)  =  100001:(100001+sz-1);   % COL1: PED (PHE.COHORTNUM?)
% ADNN(2:end , 2)  =  PHE.SRR;                % COL2: ID
% ADNN(2:end , 3)  =  zeros(sz,1);            % COL3: DAD
% ADNN(2:end , 4)  =  zeros(sz,1);            % COL4: MOM
% ADNN(2:end , 5)  =  PHE.SEX+1;              % COL5: SEX
% ADNN(2:end , 6)  =  PHE.AD+1;               % ROW1: AD (case=2;ctrl=1)


ADNN(2:end , 1)  =  PHE.SRR;                % COL1: ID
ADNN(2:end , 2)  =  PHE.AD;                 % COL2: AD
ADNN(2:end , 3)  =  PHE.COHORTNUM;          % COL3: COHORT
ADNN(2:end , 4)  =  round(PHE.AGE);         % COL4: AGE
ADNN(2:end , 5)  =  PHE.APOE;               % COL5: APOE
ADNN(2:end , 6)  =  PHE.BRAAK;              % ROW1: BRAAK

ADNN(1 , 7:end)  =  ADMX.CHRPOS';   % ROW1: CHRPOS







%------------------------------------------------------------
% CREATE  LABEL MATRIX

FULLMX = ADNN;

VX   = (FULLMX(2:end,7:end))';
TL   = (FULLMX(2:end,6))';

LX = zeros( 2 , numel(TL) );  % Make 2-D class label matrix
LX(1,:) = TL==2;              % In col-1 set CASEs index 1's
LX(2,:) = TL==1;              % In col-2 set CTRLs index 1's







%------------------------------------------------------------
% CREATE FEATURE MATRIX

ALT = zeros(size(ADMX.ALT));
ALT(ADMX.ALT == 'A') = 10;
ALT(ADMX.ALT == 'T') = 20;
ALT(ADMX.ALT == 'C') = 30;
ALT(ADMX.ALT == 'G') = 40;

REF = zeros(size(ADMX.REF));
REF(ADMX.REF == 'A') = 10;
REF(ADMX.REF == 'T') = 20;
REF(ADMX.REF == 'C') = 30;
REF(ADMX.REF == 'G') = 40;

ALTTXT = string(ALT);
ALTTXT(ALT==10) = "A";
ALTTXT(ALT==20) = "T";
ALTTXT(ALT==30) = "C";
ALTTXT(ALT==40) = "G";

REFTXT = string(REF);
REFTXT(REF==10) = "A";
REFTXT(REF==20) = "T";
REFTXT(REF==30) = "C";
REFTXT(REF==40) = "G";


UNK = REF.*0 + 1;
UNKTXT = string(REF);
UNKTXT(:) = "U";

%------------------------------------------------------------
TXTREFREF = string([zeros(size(REFTXT)); zeros(size(REFTXT))]);
TXTREFREF(1:2:end) = REFTXT;
TXTREFREF(2:2:end) = REFTXT;
TXTREFREF = TXTREFREF';
TXTREFREF = repmat(TXTREFREF,size(vMX,1),1);


TXTUNKUNK = string([zeros(size(UNKTXT)); zeros(size(UNKTXT))]);
TXTUNKUNK(1:2:end) = UNKTXT;
TXTUNKUNK(2:2:end) = UNKTXT;
TXTUNKUNK = TXTUNKUNK';
TXTUNKUNK = repmat(TXTUNKUNK,size(vMX,1),1);


TXTREFALT = string([zeros(size(REFTXT)); zeros(size(ALTTXT))]);
TXTREFALT(1:2:end) = REFTXT;
TXTREFALT(2:2:end) = ALTTXT;
TXTREFALT = TXTREFALT';
TXTREFALT = repmat(TXTREFALT,size(vMX,1),1);


TXTALTALT = string([zeros(size(ALTTXT)); zeros(size(ALTTXT))]);
TXTALTALT(1:2:end) = ALTTXT;
TXTALTALT(2:2:end) = ALTTXT;
TXTALTALT = TXTALTALT';
TXTALTALT = repmat(TXTALTALT,size(vMX,1),1);
%------------------------------------------------------------


HAPLO = [zeros(size(vMX)) zeros(size(vMX))];
HAPLO(:,1:2:end) = vMX;
HAPLO(:,2:2:end) = vMX;

HAPLOTXT = string(HAPLO);

REFREF = HAPLO == 0;
UNKUNK = HAPLO == 1;
REFALT = HAPLO == 2;
ALTALT = HAPLO == 3;


HAPLOTXT(REFREF) = TXTREFREF(REFREF);
HAPLOTXT(UNKUNK) = TXTUNKUNK(UNKUNK);
HAPLOTXT(REFALT) = TXTREFALT(REFALT);
HAPLOTXT(ALTALT) = TXTALTALT(ALTALT);


[a,b,c] = unique(HAPLOTXT);


A = strcmp(HAPLOTXT,"A");
T = strcmp(HAPLOTXT,"T");
C = strcmp(HAPLOTXT,"C");
G = strcmp(HAPLOTXT,"G");
U = strcmp(HAPLOTXT,"U");

ASUM = sum(A(:));
TSUM = sum(T(:));
CSUM = sum(C(:));
GSUM = sum(G(:));
USUM = sum(U(:));

NTOT = ASUM+TSUM+CSUM+GSUM+USUM;

APCT = ASUM / NTOT;
TPCT = TSUM / NTOT;
CPCT = CSUM / NTOT;
GPCT = GSUM / NTOT;
UPCT = USUM / NTOT;


fprintf('\n-----------------------------------\n')
fprintf('A nucleotide count: %8.f (%0.f %%)\n', ASUM, APCT*100 )
fprintf('T nucleotide count: %8.f (%0.f %%)\n', TSUM, TPCT*100 )
fprintf('C nucleotide count: %8.f (%0.f %%)\n', CSUM, CPCT*100 )
fprintf('G nucleotide count: %8.f (%0.f %%)\n', GSUM, GPCT*100 )
fprintf('U unknown bp count: %8.f (%0.f %%)\n', USUM, UPCT*100 )
fprintf('total bases  count: %8.f \n', NTOT )
fprintf('-----------------------------------\n')



% %------------------------------
% HAP = HAPLO;
% HAP(A) = 1;
% HAP(T) = 2;
% HAP(C) = 3;
% HAP(G) = 4;
% HAP(U) = 0;
% %------------------------------

%------------------------------
HAP = HAPLO;
HAP(A) = -2;
HAP(T) = -1;
HAP(C) = 1;
HAP(G) = 2;
HAP(U) = 0;
%------------------------------





NUCTXT = [string(ADNN(2:end , 1:6)) HAPLOTXT];

NUCNUM = [ADNN(2:end , 1:6) HAP];



% keyboard


%------------------------------------------------------------
% keyboard
% CP = [0 0 0 0 0 0 ADMX.CHRPOS'];
% 
% HAPLOTXT = [HAPLOTXT(1,:); HAPLOTXT];
% HAPLOTXT(1 , :) = strsplit(num2str(CP));
% 
% HAPLONUM = [HAPLONUM(1,:); HAPLONUM];
% HAPLONUM(1 , :) = CP;


% j=-1;
% k=0;
% for i = 1:(size(HAPLO,2)/2)
% j=j+2;
% k=k+2;
% 
% Hj = HAPLO(:,j);
% Hk = HAPLO(:,k);
% 
% HAPLO(Hj == 0) = REF(i);
% HAPLO(Hk == 0) = REF(i);
% HAPLO(Hj == 1) = REF(i);
% HAPLO(Hk == 1) = REF(i);
% HAPLO(Hj == 2) = REF(i);
% HAPLO(Hk == 2) = ALT(i);
% HAPLO(Hj == 3) = ALT(i);
% HAPLO(Hk == 3) = ALT(i);
% 
% 
% HAPLOTYPE(Hj == 0) = REFTXT(i);
% HAPLOTYPE(Hk == 0) = REFTXT(i);
% HAPLOTYPE(Hj == 1) = REFTXT(i);
% HAPLOTYPE(Hk == 1) = REFTXT(i);
% HAPLOTYPE(Hj == 2) = REFTXT(i);
% HAPLOTYPE(Hk == 2) = ALTTXT(i);
% HAPLOTYPE(Hj == 3) = ALTTXT(i);
% HAPLOTYPE(Hk == 3) = ALTTXT(i);
% 
% end
%------------------------------------------------------------








%------------------------------------------------------------
% ms = min(size(FULLMX));
% if ms > 9; disp(int64(FULLMX(1:7,1:7))); else disp(int64(FULLMX(1:ms,1:ms))); end
% G = sum(sum(  FULLMX(2:end , 4:end)>0, 2)>0) / (numel(TL));
% fprintf('COVERAGE: % 0.2f \n\n',G)
end
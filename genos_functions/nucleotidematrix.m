function [NUCNUM, NUCTXT] = nucleotidematrix(ADMX,CASE,CTRL,USNP,PHE,varargin)
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




% USE CUSTOM OR DEFAULT NUMERICAL VARIANT LABELS
%------------------------------------------------
if nargin > 5
    ATCGU = varargin{1};
else
           %   A  T  C  G  U
    ATCGU = [ -2  3 -4  5  0 ];
end






% CREATE VARIANT FEATURE MATRIX
%------------------------------------------------
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





% UPDATE NUMERICAL VARIANT LABELS
%------------------------------------------------
         %refref  unkunk  refalt  altalt
REFUNKALT = [-1     0       2      5];

refref = vMX==0;
unkunk = vMX==1;
refalt = vMX==2;
altalt = vMX==3;

PVX = vMX;
PVX(refref) =  REFUNKALT(1);     % HOMREF: -1
PVX(unkunk) =  REFUNKALT(2);     % UNCALL:  0
PVX(refalt) =  REFUNKALT(3);     % HETALT:  2
PVX(altalt) =  REFUNKALT(4);     % HOMALT:  5



% ADD PHENOTYPE INFO TO VARIANT MATRIX
%------------------------------------------------
PVMX = padarray(PVX,[0 9],0,'pre');
PVMX(: , 1)  =  PHE.SRR;        % COL1: ID
PVMX(: , 2)  =  PHE.AD .*2 -1;  % COL2: AD
PVMX(: , 3)  =  PHE.COHORTNUM;  % COL3: COHORT
PVMX(: , 4)  =  PHE.APOE;       % COL5: APOE
PVMX(: , 5)  =  PHE.SEX;        % COL4: SEX
PVMX(: , 6)  =  PHE.AGE;        % COL5: AGE
PVMX(: , 7)  =  PHE.BRAAK;      % COL5: BRAAK
PVMX(: , 8)  =  PHE.BRAGEZ;     % COL6: BRAGEZ
PVMX(: , 9)  =  PHE.AGEZ;       % COL4: AGEZ



% CREATE LABEL MATRIX
%------------------------------------------------
PVL = zeros( 2 , numel(PHE.AD) );  % Make 2-D class label matrix
PVL(1,:) = PHE.AD==1;              % In col-1 set CASEs index 1's
PVL(2,:) = PHE.AD==0;              % In col-2 set CTRLs index 1's

% ADNN(1 , 7:end)  =  ADMX.CHRPOS';   % ROW1: CHRPOS











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


% [a,b,c] = unique(HAPLOTXT);


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



% ASSIGN NUMERICAL VALUES TO NUCLEOTIDE LETTERS
%------------------------------------------------
HAP = HAPLO;
HAP(A) = ATCGU(1);
HAP(T) = ATCGU(2);
HAP(C) = ATCGU(3);
HAP(G) = ATCGU(4);
HAP(U) = ATCGU(5);
%------------------------------




% CREATE FINAL OUTPUT MATRICES
%------------------------------------------------
NUCTXT = [string(PVMX(: , 1:9)) HAPLOTXT];

NUCNUM = [PVMX(: , 1:9) HAP];























% DISPLAY MATRIX RANK
%------------------------------------------------

NUC = NUCNUM(:,9:end);

fprintf('\n\n---------------- MATRIX RANK --------------------\n')
fprintf('MATRIX ROWS: %4.f  (PEOPLE)\n', size(NUC,1));
fprintf('MATRIX COLS: %4.f  (VARIANTS)\n', size(NUC,2));
fprintf('MATRIX RANK: %4.f  (FULL-RANK WHEN COLS==RANK)\n\n', rank(NUC));



% DISPLAY FEATURE MATRIX PREVIEW
%------------------------------------------------
T = NUCTXT(1:9,1:15);
T = array2table(T);
T.Properties.VariableNames = {'ID','AD','COH','APOE','SEX','AGE',...
                              'BRAAK','BRAGEZ','AGEZ',...
                             ['V1a_' char(ADMX.GENE(1))],...
                             ['V1b_' char(ADMX.GENE(1))],...
                             ['V2a_' char(ADMX.GENE(2))],...
                             ['V2b_' char(ADMX.GENE(2))],...
                             ['V3a_' char(ADMX.GENE(3))],...
                             ['V3b_' char(ADMX.GENE(3))]};
T.ID     = char(T.ID);
T.AD     = str2double(T.AD);
T.COH    = str2double(T.COH);
T.APOE   = str2double(T.APOE);
T.SEX    = str2double(T.SEX);
T.AGE    = str2double(T.AGE);
T.BRAAK  = str2double(T.BRAAK);
T.BRAGEZ = str2double(T.BRAGEZ);
T.AGEZ   = str2double(T.AGEZ);
disp(T)


end
function [PVMX, VMX, DVMX, DMX, YVEC, YDUM] = mkdx(LX,CA,CO,UC,PHE,varargin)
% function [PVMX, DVMX] = mkdumx(CA,CO,UC,PHE,varargin)



SRR = PHE.SRR;

MX  = zeros( size(SRR,1) , size(CA,1) );


for nn = 1:size(CA,1)

    CASES = CA{nn};
    CTRLS = CO{nn};

    UNCAS = [UC{nn} (-ones(numel(UC{nn}),1))];
    
    CACOUN = [CASES; CTRLS; UNCAS];
    
    if any(any(CACOUN))
        CACOSRR = CACOUN(:,1);
        CACOHH  = CACOUN(:,2);
        [~,Aj] = ismember(SRR , CACOSRR );
        Af = Aj(Aj>0);
        % UNCALL:-1  HOMREF:0  HETALT:1  HOMALT:2  
        MX(Aj>0,nn) = CACOHH(Af); 
    end

end


vMX = MX + 20;   % UNCALL:19  HOMREF:20  HETALT:21  HOMALT:22
hMX = MX + 20;   % UNCALL:19  HOMREF:20  HETALT:21  HOMALT:22

% keyboard


if nargin > 5
    v = varargin{1};
else % DEFAULT
    v = [-1 0 2 4];
end




%---- MAKE vMX (ONE-COL-PER-SNP 'COOL' MX)
%----------------------------------------------
% vMX : UNKUNK:19  REFREF:20  ALTREF:21  ALTALT:22

vMX(vMX==19) = v(2);  % UNKUNK
vMX(vMX==20) = v(1);  % REFREF
vMX(vMX==21) = v(3);  % ALTREF
vMX(vMX==22) = v(4);  % ALTALT






%---- MAKE dMX (TWO-COL-PER-SNP 'HOT' MX)
%----------------------------------------------
% hMX : UNKUNK:19  REFREF:20  ALTREF:21  ALTALT:22

iMX = hMX;
jMX = hMX;

iMX(iMX==19) = v(2);  % {UNK}UNK
iMX(iMX==20) = v(1);  % {REF}REF
iMX(iMX==21) = v(3);  % {ALT}REF
iMX(iMX==22) = v(3);  % {ALT}ALT

jMX(jMX==19) = v(2);  % UNK{UNK}
jMX(jMX==20) = v(1);  % REF{REF}
jMX(jMX==21) = v(1);  % ALT{REF}
jMX(jMX==22) = v(4);  % ALT{ALT}


dMX = [hMX hMX];
dMX(:,1:2:end) = iMX;
dMX(:,2:2:end) = jMX;






%---- MAKE PVMX
%----------------------------------------------
PVMX = [zeros(size(vMX,1),9) vMX];
PVMX(: , 1)  =  PHE.SRR;        % COL1: ID
PVMX(: , 2)  =  PHE.AD;         % COL2: AD
PVMX(: , 3)  =  PHE.COHORTNUM;  % COL3: COHORT
PVMX(: , 4)  =  PHE.AGE;        % COL4: AGE
PVMX(: , 5)  =  PHE.APOE;       % COL5: APOE
PVMX(: , 6)  =  PHE.SEX;        % COL7: SEX
PVMX(: , 7)  =  PHE.BRAAK;      % COL6: BRAAK
PVMX(: , 8)  =  PHE.BRAAK;      % COL6: BRAAK
PVMX(: , 9)  =  PHE.BRAAK;      % COL6: BRAAK



DVMX = [zeros(size(dMX,1),9) dMX];
DVMX(: , 1)  =  PHE.SRR;        % COL1: ID
DVMX(: , 2)  =  PHE.AD;         % COL2: AD
DVMX(: , 3)  =  PHE.COHORTNUM;  % COL3: COHORT
DVMX(: , 4)  =  PHE.AGE;        % COL4: AGE
DVMX(: , 5)  =  PHE.APOE;       % COL5: APOE
DVMX(: , 6)  =  PHE.SEX;        % COL7: SEX
DVMX(: , 7)  =  PHE.BRAAK;      % COL6: BRAAK
DVMX(: , 8)  =  PHE.BRAAK;      % COL6: BRAAK
DVMX(: , 9)  =  PHE.BRAAK;      % COL6: BRAAK





% BRAGE = AGE - BRAD
%------------------------------------------------
% COLS:  PVMX( 1 , 2, 3 , 4 ,  5 , 6 ,  7  ,  8  ,  9    ...)
% START: PVMX(SRR,AD,COH,AGE,APOE,SEX,BRAAK,BRAAK,BRAAK.....)
% END:   PVMX(SRR,AD,COH,AGE,APOE,BRAD,BRAAK,AGEz,BRAGEz....)
% [BRAX] = bragepvmv(PVTR);

% GET BRAAK & AGE (BRAGE) WEIGHTS
PVMX = bragepvmv(PVMX);
DVMX = bragepvmv(DVMX);




% MAKE ANOTHER MATRIX THAT IS ONLY THE VARIANT COLUMNS (REMOVE PHE COLUMNS)
VMX = PVMX(:,10:end);
DMX = DVMX(:,10:end);




% CREATE AN Nx1 LABEL ARRAY
% YVEC(:,1)==1  :CASE
% YVEC(:,1)==0  :CTRL

YVEC = PVMX(:,2);


% CREATE AN Nx2 LABEL MATRIX
% YDUM(:,1)==1  :CASE
% YDUM(:,2)==1  :CTRL

YDUM = (dummyvar( (YVEC~=1)+1 ));



% keyboard
% ms = min(size(PVMX));
% if ms > 11; disp(int64(PVMX(1:11,1:11))); else disp(int64(PVMX(1:ms,1:ms))); end
% G = sum(sum(  PVMX(: , 10:end)>0, 2)>0) / (numel(TL));
% fprintf('COVERAGE: % 0.2f \n\n',G)
end
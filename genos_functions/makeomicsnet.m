function [ADNN, VX, LX] = makeomicsnet(ADMX,CASE,CTRL,USNP,PHE,varargin)


if nargin > 5
    REFUNKALT = varargin{1};
else
              %refref unkunk refalt altalt
    REFUNKALT = [-1     0     1     3];
end


CACOID = PHE.SRR;

vMX  = zeros( size(CACOID,1) , size(ADMX,1) );


for nn = 1:size(CASE,1)

    CASES = CASE{nn};
    CTRLS = CTRL{nn};

    CACO = [CASES; CTRLS];
    
    UNCAS = USNP{nn};

    if any(any(CACO) | any(UNCAS))
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
cacoMX(vMX==1) =  REFUNKALT(2);     % UNCALL: -1
cacoMX(vMX==2) =  REFUNKALT(3);     % HETALT:  1
cacoMX(vMX==3) =  REFUNKALT(4);     % HOMALT:  3

ADNN = padarray(cacoMX,[0 6],0,'pre');
ADNN(: , 1)  =  PHE.SRR;        % COL1: ID
ADNN(: , 2)  =  PHE.AD;         % COL2: AD
ADNN(: , 3)  =  PHE.COHORTNUM;  % COL3: COHORT
ADNN(: , 4)  =  PHE.AGE;        % COL4: AGE
ADNN(: , 5)  =  PHE.APOE;       % COL5: APOE
ADNN(: , 6)  =  PHE.BRAAK;      % COL6: BRAAK

VX   = (ADNN(:,7:end))';
TL   = (ADNN(:,2))';
LX = zeros( 2 , numel(TL) );  % Make 2-D class label matrix
LX(1,:) = TL==1;              % In col-1 set CASEs index 1's
LX(2,:) = TL==0;              % In col-2 set CTRLs index 1's





% ms = min(size(FULLMX));
% if ms > 9; disp(int64(FULLMX(1:7,1:7))); else disp(int64(FULLMX(1:ms,1:ms))); end
% G = sum(sum(  FULLMX(: , 7:end)>0, 2)>0) / (numel(TL));
% fprintf('COVERAGE: % 0.2f \n\n',G)
end
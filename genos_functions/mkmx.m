function [PVMX, VX, LX] = mkmx(CA,CO,UC,PHE,varargin)


SRR = PHE.SRR;

vMX  = zeros( size(SRR,1) , size(CA,1) );


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
        vMX(Aj>0,nn) = CACOHH(Af); 
    end

end


vMX = vMX + 20;   % UNCALL:19  HOMREF:20  HETALT:21  HOMALT:22



if nargin == 6
    v = varargin{1};
    vMX(vMX==19) = v(2);  % UNCALL
    vMX(vMX==20) = v(1);  % REF/REF
    vMX(vMX==21) = v(3);  % REF/ALT
    vMX(vMX==22) = v(4);  % ALT/ALT
else % DEFAULT
    vMX(vMX==19) =  0;    % UNCALL
    vMX(vMX==20) = -1;    % REF/REF
    vMX(vMX==21) =  2;    % REF/ALT
    vMX(vMX==22) =  3;    % ALT/ALT
end




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




VX   = (PVMX(:,10:end))';
TL   = (PVMX(:,2))';

LX = zeros( 2 , numel(TL) );  % Make 2-D class label matrix
LX(1,:) = TL==1;              % In col-1 set CASEs index 1's
LX(2,:) = TL==0;              % In col-2 set CTRLs index 1's


% keyboard
% ms = min(size(PVMX));
% if ms > 11; disp(int64(PVMX(1:11,1:11))); else disp(int64(PVMX(1:ms,1:ms))); end
% G = sum(sum(  PVMX(: , 10:end)>0, 2)>0) / (numel(TL));
% fprintf('COVERAGE: % 0.2f \n\n',G)
end
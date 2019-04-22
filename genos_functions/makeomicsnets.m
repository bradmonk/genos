function [FULLMX, VX, LX] = makeomicsnets(ADMX,CASE,CTRL,USNP,PHE,varargin)


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





    % Currently...
    %   HOMREF:  0
    %   UNCALL:  1
    %   HETALT:  2  
    %   HOMALT:  3
    %
    % Make...


% keyboard

    cacoMX = vMX;

    cacoMX(vMX==0) =  REFUNKALT(1);     % HOMREF: -1
    cacoMX(vMX==1) =  REFUNKALT(2);     % UNCALL: -1
    cacoMX(vMX==2) =  REFUNKALT(3);     % HETALT:  1
    cacoMX(vMX==3) =  REFUNKALT(4);     % HOMALT:  3






% ADNN = padarray(cacoMX,[1 4],0,'pre');

% ADNN(2:end , 1)  =  PHE.SRR;        % COL1: subject IDs
% ADNN(2:end , 2)  =  PHE.AD;         % COL2: AD (case=1;ctrl=0)
% ADNN(2:end , 3)  =  round(PHE.AGE); % COL3: AGE
% ADNN(2:end , 4)  =  PHE.COHORTNUM;  % COL4: COHORT
% ADNN(1 , 5:end)  =  ADMX.CHRPOS';   % ROW1: CHRPOS


% % CREATE FEATURE MATRIX AND LABEL MATRIX
% FULLMX = ADNN;
% 
% VX   = (FULLMX(2:end,5:end))';
% TL   = (FULLMX(2:end,2))';
% 
% LX = zeros( 2 , numel(TL) );  % Make 2-D class label matrix
% LX(1,:) = TL==1;              % In col-1 set CASEs index 1's
% LX(2,:) = TL==0;              % In col-2 set CTRLs index 1's
% 



ADNN = padarray(cacoMX,[0 6],0,'pre');

ADNN(: , 1)  =  PHE.SRR;        % COL1: ID
ADNN(: , 2)  =  PHE.AD;         % COL2: AD
ADNN(: , 3)  =  PHE.COHORTNUM;  % COL3: COHORT
ADNN(: , 4)  =  PHE.APOE;       % COL5: APOE
ADNN(: , 5)  =  PHE.BRAGEZ;     % COL6: BRAAK
ADNN(: , 6)  =  PHE.AGEZ;       % COL4: AGE


FULLMX = ADNN;

VX   = (FULLMX(:,7:end))';
TL   = (FULLMX(:,2))';

LX = zeros( 2 , numel(TL) );  % Make 2-D class label matrix
LX(1,:) = TL==1;              % In col-1 set CASEs index 1's
LX(2,:) = TL==0;              % In col-2 set CTRLs index 1's





ms = min(size(FULLMX));
if ms > 9; disp(int64(FULLMX(1:7,1:7))); else disp(int64(FULLMX(1:ms,1:ms))); end
G = sum(sum(  FULLMX(: , 7:end)>0, 2)>0) / (numel(TL));
fprintf('COVERAGE: % 0.2f \n\n',G)
end
function [FULLMX, VX, LX] = makennetwin_mike(ADMX,CASE,CTRL,USNP,PHE,varargin)


if nargin > 5
    doQuadratic = varargin{1};
else
    doQuadratic = 0;
end


CACOID = PHE.SRR;

vMX  = zeros( size(CACOID,1) , size(ADMX,1) );
%dims are people x snps


for nn = 1:size(CASE,1)

    CASES = CASE{nn};
    CTRLS = CTRL{nn};

    CACO = [CASES; CTRLS];
    
    UNCAS = USNP{nn};

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





if doQuadratic

    % Currently...
    %   HOMREF:  0
    %   UNCALL:  1
    %   HETALT:  2  
    %   HOMALT:  3
    %
    % Make...

    cacoMX = zeros( size(CACOID,1) , (size(ADMX,1)*3) );

    for i=1:size(vMX,1)
        for j=1:size(vMX,2)
            if vMX(i,j)==0
                cacoMX(i,(j*3)-2) = 1;
            elseif vMX(i,j)==2
                cacoMX(i,(j*3)-1) = 1;
            elseif vMX(i,j)==3
                cacoMX(i,j*3) = 1;
            end
            %cacoMX(cacoMX==0) = -2;     % HOMREF: -2
            %cacoMX(cacoMX==1) =  0;     % UNCALL:  0
            %cacoMX(cacoMX==2) =  2;     % HETALT:  2
            %cacoMX(cacoMX==3) =  60;     % HOMALT:  6
        end
    end

else

    % Currently...
    %   HOMREF:  0
    %   UNCALL:  1
    %   HETALT:  2  
    %   HOMALT:  3
    %
    % Make...
    % 
    %   HOMREF: -1
    %   UNCALL:  0
    %   HETALT:  1  
    %   HOMALT:  2

    cacoMX = vMX - 1;
end


%ADNN = padarray(cacoMX,[1 4],0,'pre');

Z = zeros(size(cacoMX,1),4);

ADNN = [Z, cacoMX];

Z = zeros(1,size(ADNN,2));

ADNN = [Z; ADNN];


ADNN(2:end , 1)  =  PHE.SRR;        % COL1: subject IDs
ADNN(2:end , 2)  =  PHE.AD;         % COL2: AD (case=1;ctrl=0)
ADNN(2:end , 3)  =  round(PHE.AGE); % COL3: AGE
ADNN(2:end , 4)  =  PHE.COHORTNUM;  % COL4: COHORT
%ADNN(1 , 5:end)  =  ADMX.CHRPOS';   % ROW1: CHRPOS






% CREATE FEATURE MATRIX AND LABEL MATRIX

FULLMX = ADNN;

VX   = (FULLMX(2:end,5:end))';
TL   = (FULLMX(2:end,2))';

LX = zeros( 2 , numel(TL) );  % Make 2-D class label matrix
LX(1,:) = TL==1;              % In col-1 set CASEs index 1's
LX(2,:) = TL==0;              % In col-2 set CTRLs index 1's










%disp(int64(FULLMX(1:10,1:10)))
G = sum(sum(  FULLMX(2:end , 4:end)>0, 2)>0) / (numel(TL));
fprintf('COVERAGE: % 0.2f \n\n',G)
end
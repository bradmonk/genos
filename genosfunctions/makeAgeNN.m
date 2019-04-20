function [ADNN, vMX] = makeAgeNN(ADMX,CASE,CTRL,USNP,PHE,varargin)


if nargin > 5
    doQuadratic = varargin{1};
else
    doQuadratic = 0;
end


CACOID = PHE.SRR;

vMX  = zeros( size(CACOID,1) , size(ADMX,1) );


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

    cacoMX = vMX;

    cacoMX(cacoMX==0) = -2;     % HOMREF: -2
    cacoMX(cacoMX==1) =  0;     % UNCALL:  0
    cacoMX(cacoMX==2) =  2;     % HETALT:  2
    cacoMX(cacoMX==3) =  6;     % HOMALT:  6

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




ADNN = padarray(cacoMX,[1 4],0,'pre');


ADNN(2:end , 1)  =  PHE.SRR;       % COL#1: subject IDs
ADNN(2:end , 2)  =  PHE.AGE;       % COL#2: AGE
ADNN(2:end , 3)  =  PHE.AD;        % COL#2: case/ctrl 1/0
ADNN(2:end , 4)  =  PHE.COHORTNUM; % COL#2: case/ctrl 1/0
ADNN(1 , 5:end)  =  ADMX.CHRPOS';  % ROW#1: CHRPOS


G = sum(sum(  ADNN(2:end , 4:end)>0, 2)>0) / (size(CACOID,1));
fprintf('COVERAGE: % 0.2f \n\n',G)



end
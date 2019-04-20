function [ADNN, caMX, coMX] = nnmake(ADMX,CASE,CTRL,USNP,PHE,varargin)


if nargin > 5
    doQuadratic = varargin{1};
else
    doQuadratic = 0;
end




CASEID = PHE.SRR(PHE.AD==1);
CTRLID = PHE.SRR(PHE.AD==0);

caMX  = zeros( size(CASEID,1) , size(ADMX,1) );
coMX  = zeros( size(CTRLID,1) , size(ADMX,1) );




% for nn = 1:size(caSNP,1)
%     caMX(:,nn) = ismember(CASEID , caSNP{nn} );
%     coMX(:,nn) = ismember(CTRLID , coSNP{nn} );
% end


for nn = 1:size(CASE,1)

%---CASES---------------------------
    CASES = CASE{nn};
    UNCAS = USNP{nn};

    if any(CASES)
        CASESRR = CASES(:,1);
        CASEHH  = CASES(:,2);

        [~,Aj] = ismember(CASEID , CASESRR );
        [~,Uj] = ismember(CASEID , UNCAS );

        Af = Aj(Aj>0);

        caMX(Aj>0,nn) = (CASEHH(Af)+1); %HETALT:2  HOMALT:3
        caMX(Uj>0,nn) = 1;              %UNCALL:1  HOMREF:0

    end
    

%---CTRLS---------------------------
    CTRLS = CTRL{nn};
    UNCON = USNP{nn};

    if any(CTRLS)
        CTRLSRR = CTRLS(:,1);
        CTRLHH  = CTRLS(:,2);

        [~,Bj] = ismember(CTRLID , CTRLSRR );
        [~,Uj] = ismember(CTRLID , UNCON );

        Bf = Bj(Bj>0);

        coMX(Bj>0,nn) = (CTRLHH(Bf) + 1);
        coMX(Uj>0,nn) = 1;

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

    cacoMX = [caMX;coMX];

    cacoMX(cacoMX==0) = -2;     % HOMREF: -2
    cacoMX(cacoMX==1) =  0;     % UNCALL:  0
    cacoMX(cacoMX==2) =  2;     % HETALT:  2
    cacoMX(cacoMX==3) =  6;     % HOMALT:  8

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

    cacoMX = [caMX;coMX] - 1;
end




ADNN = padarray(cacoMX,[1 3],0,'pre');

nc=length(CASEID)+1;


ADNN(2:end , 1)  =  [CASEID; CTRLID];   % COL#1: subject IDs
ADNN(2:nc ,  2)  =  1;                  % COL#2: case/ctrl 1/0
ADNN(2:end , 3)  =  1;                  % COL#3: bias
ADNN(1 , 4:end)  =  ADMX.CHRPOS';       % ROW#1: CHRPOS


G = sum(sum(  ADNN(2:end , 4:end)>0, 2)>0) / (size(CTRLID,1)+size(CASEID,1));
fprintf('COVERAGE: % 0.2f \n\n',G)



end
function [ADNN, caMX, coMX] = makeNNMX(ADMX,SNPCASE,SNPCTRL,CASEID,CTRLID)


[~,i] = sort(ADMX.FISHP);

ADMX     = ADMX(i,:);
SNPCASE  = SNPCASE(i);
SNPCTRL  = SNPCTRL(i);
ADMX.VID = (1:size(ADMX,1))';



VID   = ADMX.VID;
caSNP = SNPCASE(VID);
coSNP = SNPCTRL(VID);
caMX  = zeros( size(CASEID,1) , size(VID,1) );
coMX  = zeros( size(CTRLID,1) , size(VID,1) );


% for nn = 1:size(caSNP,1)
%     caMX(:,nn) = ismember(CASEID , caSNP{nn} );
%     coMX(:,nn) = ismember(CTRLID , coSNP{nn} );
% end


for nn = 1:size(caSNP,1)

%---CASES---------------------------
    CASES = caSNP{nn};

    if any(CASES)
        CASESRR = CASES(:,1);
        CASEHH  = CASES(:,2);
        [~,Aj] = ismember(CASEID , CASESRR );
        Af = Aj(Aj>0);
        caMX(Aj>0,nn) = CASEHH(Af);
    end
    

%---CTRLS---------------------------
    CTRLS = coSNP{nn};

    if any(CTRLS)
        CTRLSRR = CTRLS(:,1);
        CTRLHH = CTRLS(:,2);
        [~,Bj] = ismember(CTRLID , CTRLSRR );
        Bf = Bj(Bj>0);
        coMX(Bj>0,nn) = CTRLHH(Bf);
    end

end



cacoMX = [caMX;coMX];

ADNN = padarray(cacoMX,[1 3],0,'pre');

nc=length(CASEID)+1;


ADNN(2:end , 1)  =  [CASEID; CTRLID];   % COL#1: subject IDs
ADNN(2:nc ,  2)  =  1;                  % COL#2: case/ctrl 1/0
ADNN(2:end , 3)  =  1;                  % COL#3: bias
ADNN(1 , 4:end)  =  ADMX.VID';          % ROW#1: VID


G = sum(sum(ADNN(2:end , 4:end), 2)>0) / (size(CTRLID,1)+size(CASEID,1));
fprintf('COVERAGE: % 0.2f \n\n',G)

end
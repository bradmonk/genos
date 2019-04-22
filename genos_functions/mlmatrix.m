function ADNN = mlmatrix(CASE,CTRL,USNP,SRR,AD,COH,AGE,APOE,BRAAK)

vMX  = zeros( size(SRR,1) , size(CASE,1) );

keyboard

for nn = 1:size(CASE,1)

    CASES = CASE{nn};
    CTRLS = CTRL{nn};

    CACO = [CASES; CTRLS];
    
    UNCAS = [USNP{nn} zeros(numel(USNP{nn}),1)];
    
    CACOUN = [CASES; CTRLS; UNCAS];
    
    

    if any(any(CACO)) || any(UNCAS)
        CACOSRR = CACO(:,1);
        CACOHH  = CACO(:,2);

        [~,Aj] = ismember(SRR , CACOSRR );
        [~,Uj] = ismember(SRR , UNCAS );

        Af = Aj(Aj>0);

        vMX(Aj>0,nn) = (CACOHH(Af)+1); %HETALT:2  HOMALT:3
        vMX(Uj>0,nn) = 1;              %UNCALL:1  HOMREF:0

    end

end


vMX(vMX==0) = -1;     % HOMREF: -1
vMX(vMX==1) =  0;     % UNCALL:  0
vMX(vMX==2) =  1;     % HETALT:  1
vMX(vMX==3) =  3;     % HOMALT:  3

%ADNN = padarray(vMX,[0 6],0,'pre');
ADNN = [zeros(size(vMX,1),6) vMX];
ADNN(: , 1)  =  SRR;        % COL1: ID
ADNN(: , 2)  =  AD;         % COL2: AD
ADNN(: , 3)  =  COH;        % COL3: COHORT
ADNN(: , 4)  =  AGE;        % COL4: AGE
ADNN(: , 5)  =  APOE;       % COL5: APOE
ADNN(: , 6)  =  BRAAK;      % COL6: BRAAK

end
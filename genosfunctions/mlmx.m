function ADNN = mlmx(CA,CO,UC,SRR,AD,COH,AGE,APOE,BRAAK)

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


vMX = vMX + 1;     % UNCALL:0  HOMREF:1  HETALT:2  HOMALT:3
vMX(vMX==1) = -1;  % UNCALL:0  HOMREF:-1  HETALT:2  HOMALT:3
vMX(vMX==3) =  5;  % UNCALL:0  HOMREF:-1  HETALT:2  HOMALT:5


%ADNN = padarray(vMX,[0 6],0,'pre');
ADNN = [zeros(size(vMX,1),6) vMX];
ADNN(: , 1)  =  SRR;        % COL1: ID
ADNN(: , 2)  =  AD;         % COL2: AD
ADNN(: , 3)  =  COH;        % COL3: COHORT
ADNN(: , 4)  =  AGE;        % COL4: AGE
ADNN(: , 5)  =  APOE;       % COL5: APOE
ADNN(: , 6)  =  BRAAK;      % COL6: BRAAK

end
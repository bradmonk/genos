function [ADNN] = cohmx(NNMX,COHCASE,COHCTRL)

COHCACO = [COHCASE; COHCTRL];
COHCACO = COHCACO(:,[2 3 4 5 6 ]);
COHCACO.BRAAK(isnan(COHCACO.BRAAK)) = 1;

COHTEMP = COHCACO(1,:);
COHTEMP.SEX=0;
COHTEMP.AGE=0;
COHTEMP.APOE=0;
%COHTEMP.AUTOPSY=0;
%COHTEMP.BRAAK=0;

COHCACO = [COHTEMP; COHCACO];
COHCACO.AGE = round(COHCACO.AGE);

ADNN = [NNMX table2array(COHCACO)];

end
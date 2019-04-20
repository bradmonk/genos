function [] = fishhist(LOCI)

a = LOCI.TRFISHP;
b = LOCI.TEFISHP;
qa = quantile(a,[.001,.999]);
qb = quantile(b,[.001,.999]);

x = -log(a);
y = -log(b);
qx = quantile(x,[.001,.999]);
qy = quantile(y,[.001,.999]);

x = x(x>qx(1)&x<qx(2));
y = y(y>qy(1)&y<qy(2));

close all
fh01 = figure('Units','normalized','OuterPosition',[.03 .07 .95 .90],...
              'Color','w','MenuBar','none');
ax01 = axes('Position',[.06 .56 .4 .4],'Color','none');
ax02 = axes('Position',[.56 .56 .4 .4],'Color','none');
ax03 = axes('Position',[.06 .06 .4 .4],'Color','none');
ax04 = axes('Position',[.56 .06 .4 .4],'Color','none');

axes(ax01); histogram(a);
title('P-VALUE TRAIN TWO-TAILED')
axes(ax02); histogram(b);
title('P-VALUE TEST TWO-TAILED')

axes(ax03); histogram(x);
title('-LOGP-VALUE TRAIN TWO-TAILED')
axes(ax04); histogram(y);
title('-LOGP-VALUE TEST TWO-TAILED')
% pause(2); close all;



a = LOCI.TRHFP;
b = LOCI.TEHFP;
qa = quantile(a,[.001,.999]);
qb = quantile(b,[.001,.999]);

x = -log(a);
y = -log(b);
qx = quantile(x,[.001,.999]);
qy = quantile(y,[.001,.999]);

x = x(x>qx(1)&x<qx(2));
y = y(y>qy(1)&y<qy(2));


fh01 = figure('Units','normalized','OuterPosition',[.01 .05 .95 .90],...
              'Color','w','MenuBar','none');
ax01 = axes('Position',[.06 .56 .4 .4],'Color','none');
ax02 = axes('Position',[.56 .56 .4 .4],'Color','none');
ax03 = axes('Position',[.06 .06 .4 .4],'Color','none');
ax04 = axes('Position',[.56 .06 .4 .4],'Color','none');

axes(ax01); histogram(a);
title('P-VALUE TRAIN ONE-TAILED')
axes(ax02); histogram(b);
title('P-VALUE TEST ONE-TAILED')

axes(ax03); histogram(x);
title('-LOGP-VALUE TRAIN ONE-TAILED')
axes(ax04); histogram(y);
title('-LOGP-VALUE TEST ONE-TAILED')



%%
%------------------------------------------------------------------------
% ***  PER PERSON  ***
%-----------------------
% 
% % COMPUTE FISHERS STATISTICS FOR THE TRAINING GROUP
% [PPFISHP, PPFISHOR] = fishp_mex(LOCI.PPTRCASEREF,LOCI.PPTRCASEALT,...
%                             LOCI.PPTRCTRLREF,LOCI.PPTRCTRLALT);
% 
% LOCI.PPTRFISHP  = PPFISHP;
% LOCI.PPTRFISHOR = PPFISHOR;
% 
% 
% 
% % COMPUTE FISHERS STATISTICS FOR THE TRAINING GROUP
% [PPFISHP, PPFISHOR] = fishp_mex(LOCI.PPTECASEREF, LOCI.PPTECASEALT,...
%                             LOCI.PPTECTRLREF, LOCI.PPTECTRLALT);
% 
% LOCI.PPTEFISHP  = PPFISHP;
% LOCI.PPTEFISHOR = PPFISHOR;
%------------------------------------------------------------------------


end
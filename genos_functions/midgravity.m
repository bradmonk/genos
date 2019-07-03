function [GRAVMU] = midgravity(XAX,MISS,HIT)


MISSG = sum((XAX .* MISS));
HITG  = sum((XAX .* HIT));
MISStot = sum(MISS);
HITtot  = sum(HIT);


GRAVMU = (MISSG + HITG) ./ (MISStot + HITtot);

end
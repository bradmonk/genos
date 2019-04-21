function [] = genecards(GENE)

web(['https://www.genecards.org/cgi-bin/carddisp.pl?gene=',GENE])

end
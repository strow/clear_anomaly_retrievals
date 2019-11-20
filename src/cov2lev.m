function y = cov2lev(c,iNumLay);

m1 = (c.lev2-c.lev1)/2;
%m2 = (c.lev3-c.lev2)/2;

if nargin == 1
  pt1 = 1:97;
else
  pt1 = 1:iNumLay;
end

y1 = m1*tanh(c.width1*(pt1-c.trans1)) + (c.lev2 + c.lev1)/2;
%y2 = m2*tanh(c.width2*(pt2-c.trans2)) + (c.lev3 + c.lev2)/2;

y = y1;

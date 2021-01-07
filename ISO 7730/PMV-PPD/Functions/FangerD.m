function [j] = FangerD(in, jj)
% @ MRF
% INPUT PARAMETERS
% DESCRIPTION: FOR var = 'design'
% in is a 1-by-6 matrix which contains the input parameters
% for calculating the Air Temperature with known PMV :
% in(1,1): Metabolic rate (met)
% in(1,2): External work, normally around 0 (met)
% in(1,3): PMV, normally 0
% in(1,4): Relative humidity (0 < Rh < 1)
% in(1,5): Clothing (clo)
% in(1,6): Relative air velocity ( m/s )
% EXAMPLE
% -----------
% in =  [1.1 0 0 0.86 1 0.1];
% [TA] = FangerR(in);
% TA = 21.55
anonymf = @(j) PMVPPDd(in, j, jj);
opts = optimset('TolFun',1e-20,'Display','none');
j = fminbnd(anonymf,1,50,opts);
end
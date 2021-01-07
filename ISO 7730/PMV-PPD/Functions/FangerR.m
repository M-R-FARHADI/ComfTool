function [PMV, PPD, TA] = FangerR (in)
% @ MRF
% INPUT PARAMETERS
% DESCRIPTION: FOR var = 'rate'
% in is a 1-by-7 matrix which contains the input parameters
% for calculating the pmv and ppd:
% in(1,1): Metabolic rate (met)
% in(1,2): External work, normally around 0 (met)
% in(1,3): Radiant temperature ( C )
% in(1,4): Air Temperature ( C )
% in(1,5): Relative humidity (0 < Rh < 1)
% in(1,6): Clothing (clo)
% in(1,7): Relative air velocity ( m/s )
% EXAMPLE
% -----------
% in = [1.1 0 19.6 19.6 0.86 1 0.1];
% [PMV, PPD] = FangerR(in);
% PMV =  -0.4846 & PPD =  9.9060
psature = @ (t) 1000 .* exp(16.6536-4030.183./(t+235));
PPDD    = @ (x) 100 - 95 .* exp (-0.2179 .* x.^2 - 0.03353 .* x.^4);
mat(1:2) = in(1:2).*58;
mat(3:7) = in(3:7);
mat(6) = mat(6)*0.155;
if mat(6) <= 0.078
    fclpmv = 1 + 1.29 * mat(6);
else
    fclpmv = 1.05 + 0.645 * mat(6);
end
mat(8) = fclpmv;
opts = optimoptions('fsolve','TolFun',1e-20,'Display',...
    'none','Algorithm','trust-region-dogleg');
anonymf = @(x) pmveq(x, mat);
x = fsolve(anonymf, [30; 30], opts);
tclpmv = x(1);
hcpmv  = x(2);
q1 = (0.303 * exp(-0.036 * mat(1)) + 0.028);
PMV = q1 * (mat(1) - mat(2) - 3.05e-3 * ...
    (5733 - 6.99 * (mat(1) - mat(2)) - mat(5)...
    *psature(mat(4))) - 0.42 * (mat(1) - mat(2)...
    -58.15) - 1.7e-5 * mat(1) * (5867 - mat(5)...
    *psature(mat(4))) - 0.0014 * mat(1)*...
    (34 - mat(4)) - 3.96e-8 * fclpmv *...
    ((tclpmv + 273) ^ 4 -(mat(3) + 273) ^ 4 )...
    - fclpmv * hcpmv * (tclpmv - mat(4)));
PPD = PPDD(PMV);
TA = mat(4);
end
function [PMV] = PMVPPDd (in, j, jj)
if jj == 0
    jj = j;
end
psature = @ (t) 1000 .* exp(16.6536-4030.183./(t+235));
mat(1:2) = in(1:2).*58;
mat(3:6) = in(3:6);
mat(5) = mat(5)*0.155;
if mat(5) <= 0.078
    fclpmv = 1 + 1.29 * mat(5);
else
    fclpmv = 1.05 + 0.645 * mat(5);
end
mat(7) = fclpmv;
opts = optimoptions('fsolve','TolFun',1e-20,'Display',...
    'none','Algorithm','trust-region-dogleg');
anonymf = @(x) pmveq2(x, mat, j, jj);
x = fsolve(anonymf, [20; 20], opts);
tclpmv = x(1);
hcpmv  = x(2);
q1 = (0.303 * exp(-0.036 * mat(1)) + 0.028);
PMV = abs(q1 * (mat(1) - mat(2) - 3.05e-3 * ...
    (5733 - 6.99 * (mat(1) - mat(2)) - mat(4)...
    *psature(j)) - 0.42 * (mat(1) - mat(2)...
    -58.15) - 1.7e-5 * mat(1) * (5867 - mat(4)...
    *psature(j)) - 0.0014 * mat(1)*...
    (34 - j) - 3.96e-8 * fclpmv *...
    ((tclpmv + 273) ^ 4 -(jj + 273) ^ 4 )...
    - fclpmv * hcpmv * (tclpmv - j))) - mat(3);
end
function out = FNPS(TA)
%FNPS Summary of this function goes here
%   Detailed explanation goes here
%   FNPS calculate saturated vapor pressure KPa.

out = exp(16.6536-4030.183/(TA+235));

end


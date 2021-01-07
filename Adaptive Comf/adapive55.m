function [to, t_comf, deg_comf, ce] = adapive55(t_prevail, to, vel)
% Based on Lady-bug comfort
%     Get adaptive comfort criteria according to ASHRAE-55.
%     Note:
%         [1] ASHRAE Standard 55 (2017). Thermal Environmental Conditions
%         for Human Occupancy. Atlanta Georgoa: American Society of Heating,
%         Refrigerating and Air Conditioning Engineers.
%     Args:
%         t_prevail: The prevailing outdoor temperature [C].  For the ASHRAE-55 adaptive
%             comfort model, this is typically the average monthly outdoor temperature.
%         to: Operative temperature [C]
%         vel: Relative air velocity [m/s]
%     Returns:
%         A dictionary containing results with the following keys
%         -   to : Operative Temperature [C].
%         -   t_comf : Adaptive comfort neutral temperature (desired by occupants) [C].
%         -   deg_comf: The difference between the operative temperature (to)
%             and the adaptive comfort neutral temperature (t_comf) [C].
%             Negative values indicate cool conditions and positive values
%             indicate warm conditions.
%         -   ce -- Cooling effect as a result of elevated air speed [C]

if t_prevail < 10.
    t_prevail = 10;
elseif t_prevail > 33.5
    t_prevail = 33.5;
end
% Get the neutral temperature
neutraltemperatureashrae55 = @(x) 0.31 * x + 17.8;
% get the neutral temperature
t_comf = neutraltemperatureashrae55(t_prevail);
deg_comf = to - t_comf;
ce = 0;
if vel >= 0.6 && to >= 25
    if vel < 0.9
        ce = 1.2;
    elseif vel < 1.2
        ce = 1.8;
    elseif vel >= 1.2
        ce = 2.2;
    end
end
end
    

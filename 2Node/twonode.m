function [tempskin, tempcore, tsens, disc, SET] = twonode...
    (ta, tr, vel, rh, met, clo, wme, init_skin, init_core)
% @MRF
% J.B. Pierce two-node model of human thermoregulation.
% with some correction from ASHRAE 
% ---------------
%  _*INPUT PARAMETERS*_
% ta        : Air temperature [C]
% tr        : Mean radiant temperature [C]
% vel       : Relative air velocity [m/s]
% rh        : Relative humidity [%]
% met       : Metabolic rate [met]
% clo       : Clothing [clo]
% wme: External work [met], normally around 0 when seated
% init_skin : initial skin themprature [C]
% init_core : initial skin themprature [C]
%     _*Returns*_
% tempskin -- Skin themperature [C]
% tempcore -- Core themperature [C]
% tsens
% disc
% SET
% ---------------

% built-in functions
saturated_vapor_pressure_torr = ...
    @(t) exp(18.6686 - 4030.183 / (t + 235.0));
% initial variables.
vapor_pressure = (rh * saturated_vapor_pressure_torr(ta)) / 100;
air_velocity = max(vel, 0.1);
k_clo = 0.25;
sbc = 5.6697e-8;               % Stefan-Boltzmann constant (W/m2K4)
csw = 170;
cdil = 150;
cstr = 0.5;
bodyweight = 70;
bodysurfacearea = 1.8;
metconvertor = 58.2;
temp_skin_neutral = 33.7;      % setpoint (neutral) value for Tsk
temp_core_neutral = 36.8;      % setpoint value for Tcr
temp_body_neutral = 36.49;     % setpoint for Tb
SKBF_N = 6.3;                  % neutral value for skin_blood_flow
% INITIAL VALUES - start of 1st experiment
temp_skin = init_skin;
temp_core = init_core;
SKBF = SKBF_N;
alpha = 0.1;
esk = 0.1;
% This variable is the pressure of the atmosphere in kPa and was taken
% from the psychrometrics.js file of the CBE comfort tool.
p = 101325.0;
p = p / 1000;
pressure_in_atmospheres = p * 0.009869;
ltime = 3600;
RM = met * metconvertor;
M = met * metconvertor;
rcl = 0.155 * clo;
fcl = 1.0 + 0.15 * clo;           % body surface area due to clothing
LR = 2.2/pressure_in_atmospheres; % Lewis Relation is 2.2 at sea level
if clo <= 0
    wcrit = 0.38 * 1/(air_velocity ^ 0.29);
    icl = 1.0;
else
    wcrit = 0.59 * 1/((air_velocity)^0.08);
    icl = 0.45;
end
hc  = 3.0 * (pressure_in_atmospheres^0.53);
hcV = 8.6 * ((air_velocity * pressure_in_atmospheres)^0.53);
hc = max(hc, hcV);
% initial estimate of Tcl
hr = 4.7;
ht = hr + hc;
ra = 1.0 / (fcl * ht);       % resistance of air layer to dry heat transfer
top = (hr * tr + hc * ta) / ht;
tcl = tr + 0.1;
j = 1;
tempskin = zeros(ltime,1);
tempcore = zeros(ltime,1);

for i = 1:ltime
    j = j + 1;
    for ii = 1:10000
        % Tcl and chr are solved iteratively
        hr = 4.0 * sbc * ...
            (((tcl + tr) / 2.0 + 273.15)^ 3.0) * 0.72;
        ht = hr + hc;
        ra = 1.0 / (fcl * ht);
        top = (hr * tr + hc * ta) / ht;
        tcl1 = top + fcl * (temp_skin - top);
        if abs(tcl - tcl1) <= 0.001
            tcl = tcl1;
            break
        end
        tcl = (tcl + tcl1) / 2;
    end
    dry = (temp_skin - top) / (ra + rcl);
    hfcs = (temp_core - temp_skin) * (5.28 + 1.163 * SKBF);
    eres = 0.0023 * M * (44.0 - vapor_pressure);
    cres = 0.0014 * M * (34.0 - ta);
    scr = M - hfcs - eres - cres - wme;
    ssk = hfcs - dry - esk;
    tcsk = 0.97 * alpha * bodyweight;
    tccr = 0.97 * (1 - alpha) * bodyweight;
    dtsk = (ssk * bodysurfacearea) / (tcsk * ltime);  % deg C per ltime
    dtcr = scr * bodysurfacearea / (tccr * ltime);    % deg C per ltime
    temp_skin = temp_skin + dtsk;
    tempskin(j,1) = temp_skin;
    temp_core = temp_core + dtcr;
    tempcore(j,1) = temp_core;
    TB = alpha * temp_skin + (1 - alpha) * temp_core;
    sksig = temp_skin - temp_skin_neutral;
    if temp_skin > temp_skin_neutral
        warms = sksig;
        colds = 0;
    else
        colds = sksig;
        warms = 0;
    end
    crsig = (temp_core - temp_core_neutral);
    if temp_core > temp_core_neutral
        warmc = crsig;
        coldc = 0;
    else
        coldc =crsig;
        warmc = 0;
    end
    bdsig = TB - temp_body_neutral;
    if TB > temp_body_neutral
        warmb = bdsig;
        coldb = 0;
    else
        coldb = bdsig;
        warmb = 0;
    end
    SKBF = (SKBF_N + cdil *...
        warmc) / (1 + cstr * colds);
    if SKBF > 90.0
        SKBF = 90.0;
    end
    if SKBF < 0.5
        SKBF = 0.5;
    end
    regsw = csw * warmb * exp(warms / 10.7);
    if regsw > 500.0
        regsw = 500.0;
    end
    ersw = 0.68 * regsw;
    % evaporative resistance of air layer
    rea = 1.0 / (LR * fcl * hc);
    % evaporative resistance of clothing (icl=.45)
    recl = rcl / (LR * icl);
    emax = (saturated_vapor_pressure_torr...
        (temp_skin) - vapor_pressure) / (rea + recl);
    prsw = ersw / emax;
    pwet = 0.06 + 0.94 * prsw;
    edif = pwet * emax - ersw;
    esk = ersw + edif;
    if pwet > wcrit
        pwet = wcrit;
        prsw = wcrit / 0.94;
        ersw = prsw * emax;
        edif = 0.06 * (1.0 - prsw) * emax;
        esk = ersw + edif;
    end
    if emax < 0
        edif = 0;
        ersw = 0;
        pwet = wcrit;
        prsw = wcrit;
        esk = emax;
    end
    esk = ersw + edif;
    mshiv = 19.4 * colds * coldc;
    M = RM + mshiv;
    alpha = 0.0418 + 0.745 / (SKBF + .585);
end
% Define new heat flow terms, coeffs, and abbreviations
RN = M - wme;  % net metabolic heat production
ecomf = 0.42 * (RN - (1 * metconvertor));
if ecomf < 0.0
    ecomf = 0.0;  % from Fanger
end
emax = emax * wcrit;
% tsens is a function of TB
tbml = (.194 / 58.15) * RN + 36.301; %lower limit for evaporative regulation
tbmh = (.347 / 58.15) * RN + 36.669; %upper limit for evaporative regulation
if (TB < tbml)
    tsens = .4685 * (TB - tbml);
elseif (TB >= tbml && TB < tbmh)
    tsens = wcrit * 4.7 * (TB - tbml) / (tbmh - tbml);
elseif (TB >= tbmh)
    tsens = wcrit * 4.7 + .4685 * (TB - tbmh);
end
% disc varies with relative thermoregulatory heat strain
% Valid only when disc>0. when disc<0, disc is tsens.
disc = 4.7 * (ersw - ecomf) / (emax - ecomf - edif);
if (disc < 0)
    disc = tsens;
end
hsk = dry + esk;  % total heat loss from skin, W
W = pwet;
PSSK = exp(18.6686 - 4030.183 / (temp_skin + 235.0));
CHRS = hr;
if met < 0.85
    CHCS = 3.0;
else
    CHCS = 5.66 * (met - 0.85) ^ 0.39;
end
if CHCS < 3.0
    CHCS = 3.0;
end
CTCS = CHCS + CHRS;
RCLOS = 1.52 / ((met - wme ) + 0.6944) - 0.1835;
RCLS = 0.155 * RCLOS;
FACLS = 1.0 + k_clo * RCLOS;
FCLS = 1.0 / (1.0 + 0.155 * FACLS * CTCS * RCLOS);
IMS = 0.45;
ICLS = IMS * CHCS / CTCS * (1 - FCLS) / (CHCS / CTCS - FCLS * IMS);
RAS = 1.0 / (FACLS * CTCS);
REAS = 1.0 / (LR * FACLS * CHCS);
RECLS = RCLS / (LR * ICLS);
HD_S = 1.0 / (RAS + RCLS);
HE_S = 1.0 / (REAS + RECLS);
delta = 0.0001;
dx = 100.0;
set_old = round(temp_skin - hsk / HD_S, 2);
while abs(dx) > 0.01
    err_1 = (hsk- HD_S * (temp_skin - set_old)- ...
        W* HE_S* (PSSK - 0.5 * (exp(...
        18.6686 - 4030.183 / (set_old + 235.0)))));
    err_2 = (hsk- HD_S * (temp_skin - ...
        (set_old + delta))- W* HE_S* (PSSK - 0.5...
        * (exp(18.6686 - 4030.183 / (set_old + delta + 235.0)))));
    set = set_old - delta * err_1 / (err_2 - err_1);
    dx = set - set_old;
    set_old = set;
end
SET = set;
end
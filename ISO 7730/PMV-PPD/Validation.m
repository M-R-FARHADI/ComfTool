clc
clear
close all
%% Load datasets
valid  = load('validdata1.txt');
valid1 = load('validdata2.txt');
in = [0.016 0.12 0.2 0.3 0.4 0.6 0.8 1 1.2 1.4];
in1 = [0.001 0.2 0.3 0.5 1.5];
in3 = 12:0.05:30;
in2 = 5:0.5:45;
poolobj = parpool('local');
%% Validation 1
disp('CALCULATING ...')
tic
parfor  i = 1: size(in,2)
    input1 = [3 0 0 0.5 0.5 in(1,i)];
    [TAmet1] = FangerD(input1, 0);
    TA3(i, 1) = TAmet1;
    input2 = [2.5 0 0 0.5 0.5 in(1,i)];
    [TAmet2] = FangerD(input2, 0);
    TA25(i, 1) = TAmet2;
    input3 = [2 0 0 0.5 0.5 in(1,i)];
    [TAmet3] = FangerD(input3, 0);
    TA2(i, 1) = TAmet3;
    input4 = [1.5 0 0 0.5 0.5 in(1,i)];
    [TAmet4] = FangerD(input4, 0);
    TA15(i, 1) = TAmet4;
    input5 = [1 0 0 0.5 0.5 in(1,i)];
    [TAmet5] = FangerD(input5, 0);
    TA1(i, 1) = TAmet5;
end
toc
%% Valitation 2
disp('CALCULATING ...')
tic
parfor jj = 1 : size(in2, 2)
    input1 = [1 0 0 0.5 0.5 in1(1,1)];
    [TAmet1] = FangerD(input1, in2(1,jj));
    TA3z(jj, 1) = TAmet1;
    input2 = [1 0 0 0.5 0.5 in1(1,2)];
    [TAmet2] = FangerD(input2, in2(1,jj));
    TA25z(jj, 1) = TAmet2;
    input3 = [1 0 0 0.5 0.5 in1(1,3)];
    [TAmet3] = FangerD(input3, in2(1,jj));
    TA2z(jj, 1) = TAmet3;
    input4 = [1 0 0 0.5 0.5 in1(1,4)];
    [TAmet4] = FangerD(input4, in2(1,jj));
    TA15z(jj, 1) = TAmet4;
    input5 = [1 0 0 0.5 0.5 in1(1,5)];
    [TAmet5] = FangerD(input5, in2(1,jj));
    TA1z(jj, 1) = TAmet5;
end
toc
%% Validation 3
disp('CALCULATING ...')
tic
parfor i = 1:size(in3,2)
    input1 = [1.1 0 i i 0.86 1 0.1];
    [PMV, PPD, TA] = FangerR (input1);
    PPDv(i,1) = PPD;
    TAa(i,1)   = TA;
end
toc
%% Display
figure('units','normalized','outerposition',[0 0 1 1])
subplot(2,2,1)
hold on
plot(TA3z, in2)
plot(TA25z, in2)
plot(TA2z, in2)
plot(TA15z, in2)
plot(TA1z, in2)
plot(valid1(:,1),valid1(:,2),'*')
ylabel('MEAN RADIANT TEMP. (C)')
xlabel('AIR TEMP. (C)')
legend('V = <0.1','V = 0.2','V = 0.3','V = 0.5','V = 1.5'...
    ,'Valid.','Location','northeast')
hold off
subplot(2,2,2)
hold on
plot(TA3, in)
plot(TA25, in)
plot(TA2, in)
plot(TA15, in)
plot(TA1, in)
plot(valid(:,1),valid(:,2),'*')
xlabel('OPERATIVE TEMP. (C)')
ylabel('AIR VEL. V (m/s)')
legend('met = 3','met = 2.5','met = 2','met = 1.5','met = 1'...
    ,'Valid.','Location','northwest')
hold off
subplot(2,2,3)
plot(TAa,PPDv)
xlabel('TEMP. (C)')
ylabel('PPD')
xlim([5 40])
delete(poolobj);
clear
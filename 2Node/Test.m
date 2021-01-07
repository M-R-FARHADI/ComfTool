clc
clear
close all
%%
[tempskin, tempcore, tsens, disc, SET] = twonode...
    (25.65, 25.65, 1.7, 50, 1.5 , 0.57, 0, 40, 36.9);
AA = [tempskin, tempcore];
%% RESULT
x = 2:1:size(AA,1);
plot(x/60, AA(2:end,1), 'b')
hold on
plot(x/60, AA(2:end,2), 'r')
xlabel('Minutes of Presence')
ylabel('T(C)')
legend('Tskin','Tcore')
ylim ([33 40])
xlim ([0 60])
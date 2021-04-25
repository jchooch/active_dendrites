%% testing_area.m

clear;clc;
mus = 0:0.1:2;
tvec = 0:0.1:24000;
[currents, muvec, movave] = NoisyInput(mus,tvec);

figure(1)
tiledlayout(2,1)
nexttile
plot(tvec, currents, tvec, movave, 'r')
nexttile
plot(tvec(1:length(muvec)), muvec)

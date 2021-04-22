%Physics567 final presentation

clear;clc;
global RT RS RD CS CD Vrd Vrs theta

%% Parameters
RT = 65; % in Mohm
RS = 50; % in Mohm
RD = 43; % in Mohm
CS = 0.26; % in nF
CD = 0.12; % in nF
Vrs = -70; % in mV
Vrd = -60; % in mV
theta = -47; %in mV (soma threshold)
reset = -52; % in mV (soma reset)

% potassium current
gAHP = 4E-6 ; % in mS
EK = -90 ; % in mV
tauK = 80 ;% in ms

%calcium current
gCa = 70E-6 ;% in mS
taum = 15 ;% in ms
tauh = 80 ;% in ms

% time vector in ms
t1 = 0 ;
t2 = 20000;
dt = 0.1;
tvec = t1:dt:t2;

%vectors (vs,vd,m,h,Iahp,ICa)
vs = zeros(size(tvec));
vd = zeros(size(tvec));
vs(1) = Vrs;
vd(1) = Vrd;
m = ones(size(tvec));
h = zeros(size(tvec));
Iahp = zeros(size(tvec));
ICa = zeros(size(tvec));

Eca = 137;  % Equilibrium potential for Ca 
            % [https://www.physiologyweb.com/lecture_notes/resting_membrane_potential/resting_membrane_potential_nernst_equilibrium_potential.html]

%% Noisy input current

sigma = 0.3; %nA
tau = 3; %ms

mus = [0, 0.03, 0.06, 0.09, 0.12, 0.15,...
    0.18, 0.21, 0.24, 0.27, 0.30, 0.33, 0.36];
%mus = arange(0, 0.36, 0.03);

currents1 = [];
for mu = mus
    current = [];
    current(1) = mu;
    for i = 2:(length(tvec) / length(mus))
        current(i) = current(i-1) + ((mu - current(i-1)) / tau)*dt + ...
            sigma*normrnd(0,1)*sqrt((2*dt) / tau);
    end
    currents1 = horzcat(currents1, current);
end
for i = length(currents1):length(tvec)
    currents1(i) = mu(end);
end
movave1 = movmean(currents1, 1000);

currents2 = [];
for mu = mus
    current = [];
    current(1) = mu;
    for i = 2:(length(tvec) / length(mus))
        current(i) = current(i-1) + ((mu - current(i-1)) / tau)*dt + ...
            sigma*normrnd(0,1)*sqrt((2*dt) / tau);
    end
    currents2 = horzcat(currents2, current);
end
for i = length(currents2):length(tvec)
    currents2(i) = mu(end);
end
movave2 = movmean(currents2, 1000);


% plot current
figure(1)
tiledlayout(1,2)
nexttile
plot(tvec, currents1, 'c-', tvec, movave1, 'r--', 'LineWidth', 0.8)
nexttile
plot(tvec, currents2, 'c-', tvec, movave2, 'r--', 'LineWidth', 0.8)

Iinjs = currents1
Iinjd = currents2
            
%% iterations

for i = 1:length(tvec)
    tsum = dt*i; %tsum is current time
    if vs(i) == reset
        vd(i+1) = vd(i) + 10; %elevated by 10 mV
        Iahp(i) = gAHP*(EK - vs(i))*exp(-tsum/tauK) ; %activated when fires?
        ICa(i) = gCa*m(i)*h(i)*(Eca -vd(i)); % 10^-3V * 10^-6 A/V = 10^-9 amps
        m(i+1) = m(i) + (dt/taum)*(-m(i) + minf(vd(i)));
        h(i+1) = h(i) + (dt/tauh)*(-h(i) + hinf(vd(i)));
        vs(i+1) = vs(i) + (dt/CS)*((Vrs - vs(i))./RS + (vd(i) - vs(i))./RT + Iahp(i) + Iinjs(i));
    elseif vs(i) < theta
        Iahp(i) = gAHP*(EK - vs(i))*exp(-tsum/tauK) ; %activated when fires?
        ICa(i) = gCa*m(i)*h(i)*(Eca -vd(i));
        m(i+1) = m(i) + (dt/taum)*(-m(i) + minf(vd(i)));
        h(i+1) = h(i) + (dt/tauh)*(-h(i) + hinf(vd(i)));
        vs(i+1) = vs(i) + (dt/CS)*((Vrs - vs(i))./RS + (vd(i) - vs(i))./RT +  Iahp(i) + Iinjs(i));
        vd(i+1) = vd(i) + (dt/CD)*((Vrd - vd(i))./RD + (vs(i) - vd(i))./RT + ICa(i) + Iinjd(i));
    elseif vs(i) == 10
        vs(i+1) = reset;
        Iahp(i) = gAHP*(EK - vs(i))*exp(-tsum/tauK) ; %activated when fires?
        ICa(i) = gCa*m(i)*h(i)*(Eca -vd(i));
        m(i+1) = m(i) + (dt/taum)*(-m(i) + minf(vd(i)));
        h(i+1) = h(i) + (dt/tauh)*(-h(i) + hinf(vd(i)));
        vd(i+1) = vd(i) + (dt/CD)*((Vrd - vd(i))./RD + (vs(i) - vd(i))./RT + ICa(i) + Iinjd(i));
    else
        vs(i+1) = 10; %for 1ms
        Iahp(i) = gAHP*(EK - vs(i))*exp(-tsum/tauK) ; %activated when fires?
        ICa(i) = gCa*m(i)*h(i)*(Eca -vd(i));
        m(i+1) = m(i) + (dt/taum)*(-m(i) + minf(vd(i)));
        h(i+1) = h(i) + (dt/tauh)*(-h(i) + hinf(vd(i)));
        vd(i+1) = vd(i) + (dt/CD)*((Vrd - vd(i))./RD + (vs(i) - vd(i))./RT + ICa(i) + Iinjd(i));
    end
end

figure(2)
plot(tvec,vs(1:length(tvec)))
%xlim([0 100])
title('Somatic potential')

figure(3)
plot(tvec,vd(1:length(tvec)))
%xlim([0 100])
title('Dendritic potential')

figure(4)
plot(tvec,vs(1:length(tvec)),tvec,vd(1:length(tvec)),'r')
title('Somatic and dendritic potentials')
xlim([12000 14000])
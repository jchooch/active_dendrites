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
gAHP = 0.000004 ; % in mS
EK = -90 ; % in mV
tauK = 80 ;% in ms

%calcium current
gCa = 0.00007 ;% in mS
taum = 15 ;% in ms
tauh = 80 ;% in ms

% time vector in ms
t1 = 0 ;
t2 = 1000 ;
dt = 1;
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

Eca = Vrd;%??what's ECa

Iinjs = 2;          %nA[?]
Iinjd = 2; %stimulus in nA[?]
%% iterations

for i = 1:length(tvec)
    tsum = dt*i; %tsum is current time
    if vs(i) == reset
        vd(i+1) = vd(i) + 10; %elevated by 10 mV
        Iahp(i) = gAHP*(EK - vs(i))*exp(-tsum/tauK) ; %activated when fires?
        ICa(i) = gCa*m(i)*h(i)*(Eca -vd(i));
        m(i+1) = m(i) + (dt/taum)*(-m(i) + minf(vd(i)));
        h(i+1) = h(i) + (dt/tauh)*(-h(i) + hinf(vd(i)));
        vs(i+1) = vs(i) + (dt/CS)*((Vrs - vs(i))./RS + (vd(i) - vs(i))./RT + Iahp(i) + Iinjs);
    elseif vs(i) < theta
        Iahp(i) = gAHP*(EK - vs(i))*exp(-tsum/tauK) ; %activated when fires?
        ICa(i) = gCa*m(i)*h(i)*(Eca -vd(i));
        m(i+1) = m(i) + (dt/taum)*(-m(i) + minf(vd(i)));
        h(i+1) = h(i) + (dt/tauh)*(-h(i) + hinf(vd(i)));
        vs(i+1) = vs(i) + (dt/CS)*((Vrs - vs(i))./RS + (vd(i) - vs(i))./RT + Iahp(i) + Iinjs);
        vd(i+1) = vd(i) + (dt/CD)*((Vrd - vd(i))./RD + (vs(i) - vd(i))./RT + ICa(i) + Iinjd);
    elseif vs(i) == 10
        vs(i+1) = reset;
        Iahp(i) = gAHP*(EK - vs(i))*exp(-tsum/tauK) ; %activated when fires?
        ICa(i) = gCa*m(i)*h(i)*(Eca -vd(i));
        m(i+1) = m(i) + (dt/taum)*(-m(i) + minf(vd(i)));
        h(i+1) = h(i) + (dt/tauh)*(-h(i) + hinf(vd(i)));
        vd(i+1) = vd(i) + (dt/CD)*((Vrd - vd(i))./RD + (vs(i) - vd(i))./RT + ICa(i) + Iinjd);
    else
        vs(i+1) = 10; %for 1ms
        Iahp(i) = gAHP*(EK - vs(i))*exp(-tsum/tauK) ; %activated when fires?
        ICa(i) = gCa*m(i)*h(i)*(Eca -vd(i));
        m(i+1) = m(i) + (dt/taum)*(-m(i) + minf(vd(i)));
        h(i+1) = h(i) + (dt/tauh)*(-h(i) + hinf(vd(i)));
        vd(i+1) = vd(i) + (dt/CD)*((Vrd - vd(i))./RD + (vs(i) - vd(i))./RT + ICa(i) + Iinjd);
    end
    
end

figure(1)
plot(tvec,vs(1:length(tvec)))
xlim([0 100])
title('Somatic potential')
figure(2)
plot(tvec,vd(1:length(tvec)))
xlim([0 100])
title('Dendritic potential')


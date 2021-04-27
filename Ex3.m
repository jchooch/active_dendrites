%% Ex3.m

clear;clc;
global RT RS RD CS CD Vrd Vrs theta Eca reset

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
gAHP = 4E-3 ; % in uS
EK = -90 ; % in mV
tauK = 80 ;% in ms

%calcium current
gCa = 70E-3 ;% in uS
taum = 15 ;% in ms
tauh = 80 ;% in ms

% time vector in ms
t1 = 0 ;
t2 = 21000;
dt = 0.1;
tvec = t1:dt:t2;

%vectors (vs,vd,m,h,Iahp,ICa)
vs = zeros(size(tvec));
vd = zeros(size(tvec));
vs(1) = Vrs;
vd(1) = Vrd;
m = zeros(size(tvec));
h = zeros(size(tvec));
Iahp = zeros(size(tvec));
ICa = zeros(size(tvec));

Eca = 137;  % Equilibrium potential for Ca 
            % [https://www.physiologyweb.com/lecture_notes/resting_membrane_potential/resting_membrane_potential_nernst_equilibrium_potential.html]

%% Noisy input current

sigma = 0.3; %nA
tau = 3; %ms
mus = 0:0.1:2;
            
%% Iterations

cvs=[];
num_iters = 10;
freqs_sum = zeros(1, length(mus));
for i = 1:num_iters
    j=0;
    [currents_d, muvec_d, movave_d] = NoisyInput(mus, tvec);
    %[currents_s, muvec_s, movave_s] = NoisyInput(0 .* mus, tvec);
    %Iinjs = currents_s;
    Iinjs = zeros(1, length(currents_d));
    Iinjd = currents_d;
    
    for i = 1:length(tvec)
        if vs(i) == reset
            j = j+1;
            tsum = j*dt;
            vd(i+1) = vd(i) + 10; %elevated by 10 mV
            Iahp(i) = gAHP*(EK - vs(i))*exp(-tsum/tauK) ; %activated when fires?
            ICa(i) = gCa*m(i)*h(i)*(Eca -vd(i)); % 10^-3V * 10^-6 A/V = 10^-9 amps
            %m(i+1) = m(i) + (dt/taum)*(-m(i) + minf(vd(i)));
            %h(i+1) = h(i) + (dt/tauh)*(-h(i) + hinf(vd(i)));
            m(i+1) = m(i) + (minf(vd(i))*(dt/taum))*(1-m(i)) - ((1 - minf(vd(i)))*(dt/taum))*m(i);
            h(i+1) = h(i) + (hinf(vd(i))*(dt/tauh))*(1-h(i)) - ((1 - hinf(vd(i)))*(dt/tauh))*h(i);
            vs(i+1) = vs(i) + (dt/CS)*((Vrs - vs(i))./RS + (vd(i) - vs(i))./RT + Iahp(i) + Iinjs(i));  
        elseif vs(i) < theta && j==0
            tsum=0;
            Iahp(i) = 0; %gAHP*(EK - vs(i))*exp(-tsum/tauK) ; %activated when fires?
            ICa(i) = gCa*m(i)*h(i)*(Eca -vd(i));        % nA
            %m(i+1) = m(i) + (dt/taum)*(-m(i) + minf(vd(i)));
            %h(i+1) = h(i) + (dt/tauh)*(-h(i) + hinf(vd(i)));
            m(i+1) = m(i) + (minf(vd(i))*(dt/taum))*(1-m(i)) - ((1 - minf(vd(i)))*(dt/taum))*m(i);
            h(i+1) = h(i) + (hinf(vd(i))*(dt/tauh))*(1-h(i)) - ((1 - hinf(vd(i)))*(dt/tauh))*h(i);
            vs(i+1) = vs(i) + (dt/CS)*((Vrs  - vs(i))./RS + (vd(i) - vs(i))./RT +  Iahp(i) + Iinjs(i));
            vd(i+1) = vd(i) + (dt/CD)*((Vrd - vd(i))./RD + (vs(i) - vd(i))./RT + ICa(i) + Iinjd(i));    
       elseif vs(i) < theta && j>0
            j =j+1;
            tsum = j*dt;
            Iahp(i) =gAHP*(EK - vs(i))*exp(-tsum/tauK) ; 
            ICa(i) = gCa*m(i)*h(i)*(Eca -vd(i));
            %m(i+1) = m(i) + (dt/taum)*(-m(i) + minf(vd(i)));
            %h(i+1) = h(i) + (dt/tauh)*(-h(i) + hinf(vd(i)));
            m(i+1) = m(i) + (minf(vd(i))*(dt/taum))*(1-m(i)) - ((1 - minf(vd(i)))*(dt/taum))*m(i);
            h(i+1) = h(i) + (hinf(vd(i))*(dt/tauh))*(1-h(i)) - ((1 - hinf(vd(i)))*(dt/tauh))*h(i);
            vs(i+1) = vs(i) + (dt/CS)*((Vrs - vs(i))./RS + (vd(i) - vs(i))./RT +  Iahp(i) + Iinjs(i));
            vd(i+1) = vd(i) + (dt/CD)*((Vrd - vd(i))./RD + (vs(i) - vd(i))./RT + ICa(i) + Iinjd(i));  
        elseif vs(i) == 10
            j = j+1;
            tsum = j*dt;
            vs(i+1) = reset;
            Iahp(i) = gAHP*(EK - vs(i))*exp(-tsum/tauK) ; %activated when fires?
            ICa(i) = gCa*m(i)*h(i)*(Eca -vd(i));
            %m(i+1) = m(i) + (dt/taum)*(-m(i) + minf(vd(i)));
            %h(i+1) = h(i) + (dt/tauh)*(-h(i) + hinf(vd(i)));
            m(i+1) = m(i) + (minf(vd(i))*(dt/taum))*(1-m(i)) - ((1 - minf(vd(i)))*(dt/taum))*m(i);
            h(i+1) = h(i) + (hinf(vd(i))*(dt/tauh))*(1-h(i)) - ((1 - hinf(vd(i)))*(dt/tauh))*h(i);
            vd(i+1) = vd(i) + (dt/CD)*((Vrd - vd(i))./RD + (vs(i) - vd(i))./RT + ICa(i) + Iinjd(i));
        else
            j=1;
            tsum = j*dt;
            vs(i+1) = 10; %for 1ms
            Iahp(i) = gAHP*(EK - vs(i))*exp(-tsum/tauK) ; %activated when fires?
            ICa(i) = gCa*m(i)*h(i)*(Eca -vd(i));
            %m(i+1) = m(i) + (dt/taum)*(-m(i) + minf(vd(i)));
            %h(i+1) = h(i) + (dt/tauh)*(-h(i) + hinf(vd(i)));
            m(i+1) = m(i) + (minf(vd(i))*(dt/taum))*(1-m(i)) - ((1 - minf(vd(i)))*(dt/taum))*m(i);
            h(i+1) = h(i) + (hinf(vd(i))*(dt/tauh))*(1-h(i)) - ((1 - hinf(vd(i)))*(dt/tauh))*h(i);
            vd(i+1) = vd(i) + (dt/CD)*((Vrd - vd(i))./RD + (vs(i) - vd(i))./RT + ICa(i) + Iinjd(i));   
        end
    end
    window_size = floor(length(tvec) / length(mus));
    num_spikes = [];
    cv = [];
    loc = cell(size(mus));
    isi = cell(size(mus));
    for i = 0:(length(mus)-1)
        num_spikes(i+1) = nnz(vs((i*window_size+1):((i+1)*window_size+1)) == 10);
        %num_spikes(i+1) = length(findpeaks(vs((i*window_size+1):((i+1)*window_size+1))))
        loc{i+1} = find(vs((i*window_size+1):((i+1)*window_size+1)) == 10);
        isi{i+1} = diff(loc{i+1}).*dt;% in ms
        cv(i+1) = std(isi{i+1})/mean(isi{i+1}); %getCV(isi{i+1});
    end
    freqs = (num_spikes ./ (window_size/1000));
    cvs = vertcat(cvs,cv);
    freqs_sum = freqs_sum + freqs;
end
ICa_ca = ICa;
ave_freqs_ca = freqs_sum ./ num_iters;
cvs_for_plotting_ca = nanmean(cvs, 1);
cvs_for_plotting_ca(isnan(cvs_for_plotting_ca)) = 0;
vs_ca = vs;
vd_ca = vd;
m_ca = m;
h_ca = h;
Iahp_ca = Iahp;

cvs=[];
freqs_sum = zeros(1, length(mus));
gCa = 0;
vs(1) = Vrs;
vd(1) = Vrd;
for i = 1:num_iters
    j=0;
    [currents_d, muvec_d, movave_d] = NoisyInput(mus, tvec);
    %[currents_s, muvec_s, movave_s] = NoisyInput(0 .* mus, tvec);
    %Iinjs = currents_s;
    Iinjs = zeros(1, length(currents_d));
    Iinjd = currents_d;
    
    for i = 1:length(tvec)
        if vs(i) == reset
            j = j+1;
            tsum = j*dt;
            Iahp(i) = gAHP*(EK - vs(i))*exp(-tsum/tauK) ; %activated when fires?
            ICa(i) = gCa*m(i)*h(i)*(Eca -vd(i)); % 10^-3V * 10^-6 A/V = 10^-9 amps
            %m(i+1) = m(i) + (dt/taum)*(-m(i) + minf(vd(i)));
            %h(i+1) = h(i) + (dt/tauh)*(-h(i) + hinf(vd(i)));
            m(i+1) = m(i) + (minf(vd(i))*(dt/taum))*(1-m(i)) - ((1 - minf(vd(i)))*(dt/taum))*m(i);
            h(i+1) = h(i) + (hinf(vd(i))*(dt/tauh))*(1-h(i)) - ((1 - hinf(vd(i)))*(dt/tauh))*h(i);
            vs(i+1) = vs(i) + (dt/CS)*((Vrs - vs(i))./RS + (vd(i) - vs(i))./RT + Iahp(i) + Iinjs(i));
            vd(i+1) = vd(i) + 10; %elevated by 10 mV
        elseif vs(i) < theta && j==0
            tsum=0;
            Iahp(i) = 0; %gAHP*(EK - vs(i))*exp(-tsum/tauK) ; %activated when fires?
            ICa(i) = gCa*m(i)*h(i)*(Eca -vd(i));
            %m(i+1) = m(i) + (dt/taum)*(-m(i) + minf(vd(i)));
            %h(i+1) = h(i) + (dt/tauh)*(-h(i) + hinf(vd(i)));
            m(i+1) = m(i) + (minf(vd(i))*(dt/taum))*(1-m(i)) - ((1 - minf(vd(i)))*(dt/taum))*m(i);
            h(i+1) = h(i) + (hinf(vd(i))*(dt/tauh))*(1-h(i)) - ((1 - hinf(vd(i)))*(dt/tauh))*h(i);
            vs(i+1) = vs(i) + (dt/CS)*((Vrs - vs(i))./RS + (vd(i) - vs(i))./RT +  Iahp(i) + Iinjs(i));
            vd(i+1) = vd(i) + (dt/CD)*((Vrd - vd(i))./RD + (vs(i) - vd(i))./RT + ICa(i) + Iinjd(i));    
       elseif vs(i) < theta && j>0
            j =j+1;
            tsum = j*dt;
            Iahp(i) =gAHP*(EK - vs(i))*exp(-tsum/tauK) ; 
            ICa(i) = gCa*m(i)*h(i)*(Eca -vd(i));
            %m(i+1) = m(i) + (dt/taum)*(-m(i) + minf(vd(i)));
            %h(i+1) = h(i) + (dt/tauh)*(-h(i) + hinf(vd(i)));
            m(i+1) = m(i) + (minf(vd(i))*(dt/taum))*(1-m(i)) - ((1 - minf(vd(i)))*(dt/taum))*m(i);
            h(i+1) = h(i) + (hinf(vd(i))*(dt/tauh))*(1-h(i)) - ((1 - hinf(vd(i)))*(dt/tauh))*h(i);
            vs(i+1) = vs(i) + (dt/CS)*((Vrs - vs(i))./RS + (vd(i) - vs(i))./RT +  Iahp(i) + Iinjs(i));
            vd(i+1) = vd(i) + (dt/CD)*((Vrd - vd(i))./RD + (vs(i) - vd(i))./RT + ICa(i) + Iinjd(i));  
        elseif vs(i) == 10
            j = j+1;
            tsum = j*dt;
            Iahp(i) = gAHP*(EK - vs(i))*exp(-tsum/tauK) ; %activated when fires?
            ICa(i) = gCa*m(i)*h(i)*(Eca -vd(i));
            %m(i+1) = m(i) + (dt/taum)*(-m(i) + minf(vd(i)));
            %h(i+1) = h(i) + (dt/tauh)*(-h(i) + hinf(vd(i)));
            m(i+1) = m(i) + (minf(vd(i))*(dt/taum))*(1-m(i)) - ((1 - minf(vd(i)))*(dt/taum))*m(i);
            h(i+1) = h(i) + (hinf(vd(i))*(dt/tauh))*(1-h(i)) - ((1 - hinf(vd(i)))*(dt/tauh))*h(i);
            vs(i+1) = reset;
            vd(i+1) = vd(i) + (dt/CD)*((Vrd - vd(i))./RD + (vs(i) - vd(i))./RT + ICa(i) + Iinjd(i));
        else
            j=1;
            tsum = j*dt;
            Iahp(i) = gAHP*(EK - vs(i))*exp(-tsum/tauK) ; %activated when fires?
            ICa(i) = gCa*m(i)*h(i)*(Eca -vd(i));
            %m(i+1) = m(i) + (dt/taum)*(-m(i) + minf(vd(i)));
            %h(i+1) = h(i) + (dt/tauh)*(-h(i) + hinf(vd(i)));
            m(i+1) = m(i) + (minf(vd(i))*(dt/taum))*(1-m(i)) - ((1 - minf(vd(i)))*(dt/taum))*m(i);
            h(i+1) = h(i) + (hinf(vd(i))*(dt/tauh))*(1-h(i)) - ((1 - hinf(vd(i)))*(dt/tauh))*h(i);
            vs(i+1) = 10; %for 1ms
            vd(i+1) = vd(i) + (dt/CD)*((Vrd - vd(i))./RD + (vs(i) - vd(i))./RT + ICa(i) + Iinjd(i));   
        end
    end
    window_size = floor(length(tvec) / length(mus));
    num_spikes = [];
    cv = [];
    loc = cell(size(mus));
    isi = cell(size(mus));
    for i = 0:(length(mus)-1)
        num_spikes(i+1) = nnz(vs((i*window_size+1):((i+1)*window_size+1)) == 10);
        %num_spikes(i+1) = length(findpeaks(vs((i*window_size+1):((i+1)*window_size+1))))
        loc{i+1} = find(vs((i*window_size+1):((i+1)*window_size+1)) == 10);
        isi{i+1} = diff(loc{i+1}).*dt;% in ms
        cv(i+1) = std(isi{i+1})/mean(isi{i+1}); %getCV(isi{i+1});
    end
    freqs = (num_spikes ./ (window_size/1000));
    cvs = vertcat(cvs,cv);
    freqs_sum = freqs_sum + freqs;
end
ICa_cd = ICa;
ave_freqs_cd = freqs_sum ./ num_iters;
cvs_for_plotting_cd = nanmean(cvs, 1);
cvs_for_plotting_cd(isnan(cvs_for_plotting_cd)) = 0;
vs_cd = vs;
vd_cd = vd;
m_cd = m;
h_cd = h;
Iahp_cd = Iahp;

%% Plots

%{
Figures, in order of production below:
From paper:
- Fig 3a fI curves for Cd2+ and Ca2+
- Fig 3b CVs of ISIs
- Fig 3c Somatic potential with Ca2+ in zoomed window
- Fig 3d Somatic potential with Cd2+ in zoomed window
- Fig 3e Zoom on 3c
- Fig 3f Gain reduction HOW???
Supplementary (ours):
- 
%}

times_for_plotting = tvec ./ 10; % times in ms for plotting

% Fig 3A ? fI curves for Cd2+ and Ca2+
figure(1)
rectangle('Position',[1200 0 200 18], 'FaceColor', [0.9 0.9 0.9])
hold on
rectangle('Position',[1400 0 700 18], 'FaceColor', [1 0.9 0.9])
hold on
plot(mus .* 1000, ave_freqs_ca, 'ko-', 'MarkerFaceColor', 'k', 'MarkerSize', 10)
hold on
plot(mus .* 1000, ave_freqs_cd, 'k--')
hold off
title('Dendritic fI curves with and without calcium')
xlim([0 2100])
legend('Ca2+', 'Cd2+')
xlabel('Mean current (pA)')
ylabel('Mean spike rate (AP/s)')

figure(2)
plot(tvec, ICa_ca .* 1000, 'k-', tvec, ICa_cd .* 1000, 'r-')
title('Calcium current over time')
xlabel('Time (ms)')
ylabel('Current (pA)')

figure(3)
tiledlayout(2,1)
nexttile
plot(tvec, vs_ca(1:length(tvec)), tvec, vd_ca(1:length(tvec)))
title('Without calcium channel blocker')
legend('v_s', 'v_d')
ylabel('Potential (mV)')
ylim([-80 40])
nexttile
plot(tvec, vs_cd(1:length(tvec)), tvec, vd_cd(1:length(tvec)))
title('With calcium channel blocker')
legend('v_s', 'v_d')
ylabel('Potential (mV)')
ylim([-80 40])
xlabel('Time (ms)')

figure(4)
tiledlayout(2,2)
nexttile
plot(tvec, m_ca(1:length(tvec)), '-', 'Color', [.5 0 .5])
title('m gating variable over time')
nexttile
plot(tvec, h_ca(1:length(tvec)), 'g-')
title('h gating variable over time')
nexttile
plot(tvec, m_cd(1:length(tvec)), '-', 'Color', [.5 0 .5])
xlabel('Time (ms)')
nexttile
plot(tvec, h_cd(1:length(tvec)), 'g-')
xlabel('Time (ms)')

figure(5)
tiledlayout(2,1)
nexttile
plot(tvec, Iahp_ca)
nexttile
plot(tvec, Iahp_cd)

figure(6)
rectangle('Position',[1200 0 200 4.5], 'FaceColor', [0.9 0.9 0.9])
hold on
rectangle('Position',[1400 0 700 4.5], 'FaceColor', [1 0.9 0.9])
hold on
plot(mus .* 1000, cvs_for_plotting_ca, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 10)
hold on
plot(mus .* 1000, cvs_for_plotting_cd, 'ko', 'MarkerSize', 10)
hold off
xlim([0 2100])
title('CVs of ISIs for different mean dendritic input currents')
xlabel('Mean current (pA)')
ylabel('CV of ISIs')

figure(7)
tiledlayout(1,2)
nexttile
plot(tvec(120001:140001), vs_ca(120001:140001))
title('Without calcium channel blocker')
xlabel('Time')
ylabel('Somatic potential (mV)')
nexttile
plot(tvec(120001:140001), vs_cd(120001:140001))
title('Without calcium channel blocker')
xlabel('Time (ms)')

%{
% Fig 1B - Stimulate soma only
figure(1)
tiledlayout(2,1)
nexttile
plot(tvec,vs(1:length(tvec)),'k-')
ylabel('Somatic voltage (mV)')
title('Fig 1B: Somatic current injection')
nexttile
plot(tvec,Iinjs,'k-', tvec(1:length(muvec_s)), muvec_s,'w-','LineWidth', 1)
%plot(tvec,Iinjs,'k-', tvec, movave1,'b--','LineWidth', 0.8)
ylabel('Noisy somatic input (pA)')
xlabel('Time (ms)')
%}


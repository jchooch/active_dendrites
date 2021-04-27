%Physics567 final presentation fig4 - figure4. run for s+d

clear;clc;
global RT RS RD CS CD Vrd Vrs theta Eca reset

%% Parameters
RT = 65; % in Mohm %model distal dendrite (original 65)
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
t2 = 48000; 
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

% mus = [0, 0.03, 0.06, 0.09, 0.12, 0.15,...
%     0.18, 0.21, 0.24, 0.27, 0.30, 0.33]; %in nA
% test soma threshold
mus_s = ones(1,12).*0.25;
mus = 0:0.1:1.1;%step 100pA
% injd
% mus = 1:0.06:2.3;
            
%% iterations
cvs=[];
num_iters = 1;
freqs_sum = zeros(1, length(mus));
for i = 1:num_iters
    j=0;
    %generate noisy current imput
    currents1 = [];
    muvec = [];
    for mu = mus_s
        current = [];
        current(1) = mu;
        for i = 2:(length(tvec) / length(mus))
            current(i) = current(i-1) + ((mu - current(i-1)) / tau)*dt + ...
                sigma*normrnd(0,1)*sqrt((2*dt) / tau);
           
        end
        muv = ones(size(current)).*mu;
        currents1 = horzcat(currents1, current);
        muvec = horzcat(muvec,muv);
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
    
    Iinjs = currents1;
    Iinjd = currents2;
    %Iinjs = zeros(1, length(currents1));
    
    for i = 1:length(tvec)
        if vs(i) == reset
            j = j+1;
            tsum = j*dt;
            vd(i+1) = vd(i) + 10; %elevated by 10 mV
            Iahp(i) = gAHP*(EK - vs(i))*exp(-tsum/tauK) ; %activated when fires?
            ICa(i) = gCa*m(i)*h(i)*(Eca -vd(i)); % 10^-3V * 10^-6 A/V = 10^-9 amps
            m(i+1) = m(i) + (dt/taum)*(-m(i) + minf(vd(i)));
            h(i+1) = h(i) + (dt/tauh)*(-h(i) + hinf(vd(i)));
            vs(i+1) = vs(i) + (dt/CS)*((Vrs - vs(i))./RS + (vd(i) - vs(i))./RT + Iahp(i) + Iinjs(i));  
        elseif vs(i) < theta && j==0
            tsum=0;
            Iahp(i) = 0; %gAHP*(EK - vs(i))*exp(-tsum/tauK) ; %activated when fires?
            ICa(i) = gCa*m(i)*h(i)*(Eca -vd(i));
            m(i+1) = m(i) + (dt/taum)*(-m(i) + minf(vd(i)));
            h(i+1) = h(i) + (dt/tauh)*(-h(i) + hinf(vd(i)));
            vs(i+1) = vs(i) + (dt/CS)*((Vrs - vs(i))./RS + (vd(i) - vs(i))./RT +  Iahp(i) + Iinjs(i));
            vd(i+1) = vd(i) + (dt/CD)*((Vrd - vd(i))./RD + (vs(i) - vd(i))./RT + ICa(i) + Iinjd(i));    
       elseif vs(i) < theta && j>0
            j =j+1;
            tsum = j*dt;
            Iahp(i) =gAHP*(EK - vs(i))*exp(-tsum/tauK) ; 
            ICa(i) = gCa*m(i)*h(i)*(Eca -vd(i));
            m(i+1) = m(i) + (dt/taum)*(-m(i) + minf(vd(i)));
            h(i+1) = h(i) + (dt/tauh)*(-h(i) + hinf(vd(i)));
            vs(i+1) = vs(i) + (dt/CS)*((Vrs - vs(i))./RS + (vd(i) - vs(i))./RT +  Iahp(i) + Iinjs(i));
            vd(i+1) = vd(i) + (dt/CD)*((Vrd - vd(i))./RD + (vs(i) - vd(i))./RT + ICa(i) + Iinjd(i));  
        elseif vs(i) == 10
            j = j+1;
            tsum = j*dt;
            vs(i+1) = reset;
            Iahp(i) = gAHP*(EK - vs(i))*exp(-tsum/tauK) ; %activated when fires?
            ICa(i) = gCa*m(i)*h(i)*(Eca -vd(i));
            m(i+1) = m(i) + (dt/taum)*(-m(i) + minf(vd(i)));
            h(i+1) = h(i) + (dt/tauh)*(-h(i) + hinf(vd(i)));
            vd(i+1) = vd(i) + (dt/CD)*((Vrd - vd(i))./RD + (vs(i) - vd(i))./RT + ICa(i) + Iinjd(i));
        else
            j=1;
            tsum = j*dt;
            vs(i+1) = 10; %for 1ms
            Iahp(i) = gAHP*(EK - vs(i))*exp(-tsum/tauK) ; %activated when fires?
            ICa(i) = gCa*m(i)*h(i)*(Eca -vd(i));
            m(i+1) = m(i) + (dt/taum)*(-m(i) + minf(vd(i)));
            h(i+1) = h(i) + (dt/tauh)*(-h(i) + hinf(vd(i)));
            vd(i+1) = vd(i) + (dt/CD)*((Vrd - vd(i))./RD + (vs(i) - vd(i))./RT + ICa(i) + Iinjd(i));   
        end
    end
    window_size = floor(length(tvec) / length(mus)); %step duration of 2s
    num_spikes_sd = []; %soma

    cv = [];
    loc_s = cell(size(mus));
    isi = cell(size(mus));
    for i = 0:length(mus)-1
        num_spikes_sd(i+1) = nnz(vs((i*window_size+1):((i+1)*window_size+1)) == 10);
 
        %num_spikes(i+1) = length(findpeaks(vs((i*window_size+1):((i+1)*window_size+1))))
        loc_s{i+1} = find(vs((i*window_size+1):((i+1)*window_size+1)) == 10);
        isi{i+1} = diff(loc_s{i+1}).*dt;% in ms
        cv(i+1) = std(isi{i+1})/mean(isi{i+1});%getCV(isi{i+1});
    end
    freqs_sd = ((num_spikes_sd * 1000) ./ (window_size*dt));
    cvs = vertcat(cvs,cv);
    freqs_sum = freqs_sum + freqs_sd;
end
ave_freqs_sd = freqs_sum ./ num_iters; %soma 
mus_sd = mus;

%save voltages 
vs_sd = vs;
vd_sd = vd;
cvs_for_plotting_sd = nanmean(cvs, 1); %soma_cv
cvs_for_plotting_sd(isnan(cvs_for_plotting_sd)) = 0;
figure(8) %%CV of isi 
plot(mus_sd,cvs_for_plotting_sd,'ok','MarkerFaceColor','r')

%% plot spike trains
figure(1)
plot(tvec(1:4*10^5),vs_sd(1:4*10^5),'k-')
hold on
plot(tvec(1:2*10^4),vs_sd(1:2*10^4),'b-') %4D
plot(tvec(15*10^4:17*10^4),vs_sd(15*10^4:17*10^4),'b-') %4E
plot(tvec(3.8*10^5:4*10^5),vs_sd(3.8*10^5:4*10^5),'b-')
plot([-2000; -2000], [-40; 0], '-k',  [10000; 20000], [-75; -75], '-k', 'LineWidth', 2)
text(-2600,-20, '40 mV', 'HorizontalAlignment','right')
text(15000,-78, '10 s', 'HorizontalAlignment','center')
set(gca, 'Visible', 'off')
legend('Somatic spike train','Location','Southeast')
hold off


%fig 4F final 2s 
figure(2) 
plot(tvec(3.8*10^5:4*10^5),vs_sd(3.8*10^5:4*10^5),'k-') %F
hold on
plot(tvec(386300:386300+1500),vs_sd(386300:386300+1500),'b-') %4F1
plot(tvec(3.8*10^5:4*10^5),vd_sd(3.8*10^5:4*10^5),'r-')
% plot([1200; 1700], [-10; -10], '-k',  [1300; 1300], [-15; -35], '-k', 'LineWidth', 2)
% text(1480,-25, '20 mV', 'HorizontalAlignment','right')
% text(1460,-7, '500 ms', 'HorizontalAlignment','center')
set(gca, 'Visible', 'off')
legend('Continued somatic spiking','Location','Southeast')
hold off

%E1 Done
figure(3)
plot(tvec(386300:386300+1500),vs_sd(386300:386300+1500),'k-') %4F1
hold on
plot(tvec(386300:386300+1500),vd_sd(386300:386300+1500),'r-')
% plot([220, 270], [-10; -10], '-k',  [220,220], [-15; -35], '-k', 'LineWidth', 2)
% text(235,-25, '20 mV', 'HorizontalAlignment','right')
% text(245,-7, '50 ms', 'HorizontalAlignment','center')
set(gca, 'Visible', 'off')
legend('Continued somatic spiking','Location','Southeast')
hold off





% %fig 4D first 2s Done
% 
% figure(2) 
% plot(tvec(1:2*10^4),vs_sd(1:2*10^4),'k-')
% hold on
% plot(tvec(1500:3000),vs_sd(1500:3000),'b-','Linewidth',1) %4D1
% plot(tvec(1:2*10^4),vd_sd(1:2*10^4),'r-')
% plot([1200; 1700], [-10; -10], '-k',  [1300; 1300], [-15; -35], '-k', 'LineWidth', 2)
% text(1480,-25, '20 mV', 'HorizontalAlignment','right')
% text(1460,-7, '500 ms', 'HorizontalAlignment','center')
% set(gca, 'Visible', 'off')
% legend('Single somatic APs','Location','Northeast')
% hold off
% 
%D1 Done
% figure()
% plot(tvec(1500:3000),vs_sd(1500:3000),'k-') %4D1
% hold on
% plot(tvec(1500:3000),vd_sd(1500:3000),'r-')
% plot([220, 270], [-10; -10], '-k',  [220,220], [-15; -35], '-k', 'LineWidth', 2)
% text(235,-25, '20 mV', 'HorizontalAlignment','right')
% text(245,-7, '50 ms', 'HorizontalAlignment','center')
% set(gca, 'Visible', 'off')
% legend('Single somatic APs','Location','Northeast')
% hold off

% %fig 4E mid 2s Done
% figure(2) 
% plot(tvec(15*10^4:17*10^4),vs_sd(15*10^4:17*10^4),'k-') %E
% hold on
% plot(tvec(152400:152400+1500),vs_sd(152400:152400+1500),'b-') %4E1
% plot(tvec(15*10^4:17*10^4),vd_sd(15*10^4:17*10^4),'r-')
% % plot([1200; 1700], [-10; -10], '-k',  [1300; 1300], [-15; -35], '-k', 'LineWidth', 2)
% % text(1480,-25, '20 mV', 'HorizontalAlignment','right')
% % text(1460,-7, '500 ms', 'HorizontalAlignment','center')
% set(gca, 'Visible', 'off')
% legend('Somatic bursting','Location','Southeast')
% hold off
% 
% %E1 Done
% figure(3)
% plot(tvec(152400:152400+1500),vs_sd(152400:152400+1500),'k-') %4E1
% hold on
% plot(tvec(152400:152400+1500),vd_sd(152400:152400+1500),'r-')
% % plot([220, 270], [-10; -10], '-k',  [220,220], [-15; -35], '-k', 'LineWidth', 2)
% % text(235,-25, '20 mV', 'HorizontalAlignment','right')
% % text(245,-7, '50 ms', 'HorizontalAlignment','center')
% set(gca, 'Visible', 'off')
% legend('Somatic bursting','Location','Northeast')
% hold off


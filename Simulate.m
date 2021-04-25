%% Simulate.m

function [ave_freqs, cvs_for_plotting] = Simulate(stim_s, stim_d, mus_soma, mus_dend, num_iters)

    if stim_s == True
        mus_s = mus_soma
    else
        mus_s = 0
    end
    cvs=[];
    freqs_sum = zeros(1, length(mus));
    for i = 1:num_iters
        j=0;
        [currents_s, muvec_s, movave_s] = NoisyInput(mus, tvec);
        [currents_d, muvec_d, movave_d] = NoisyInput(mus, tvec);
        Iinjs = currents_s;
        Iinjd = zeros(1, length(currents_d))

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
        window_size = floor(length(tvec) / length(mus));
        num_spikes = [];
        cv = [];
        loc = cell(size(mus));
        isi = cell(size(mus));
        for i = 0:0:(length(mus)-1)
            num_spikes(i+1) = nnz(vs((i*window_size+1):((i+1)*window_size+1)) == 10);
            %num_spikes(i+1) = length(findpeaks(vs((i*window_size+1):((i+1)*window_size+1))))
            loc{i+1} = find(vs((i*window_size+1):((i+1)*window_size+1)) == 10);
            isi{i+1} = diff(loc{i+1}).*dt;% in ms
            cv(i+1) = std(isi{i+1})/mean(isi{i+1}); %getCV(isi{i+1});
        end
        freqs = ((num_spikes * 1000) ./ (window_size*dt));
        cvs = vertcat(cvs,cv);
        freqs_sum = freqs_sum + freqs;
    end
    
    ave_freqs = freqs_sum ./ num_iters;
    cvs_for_plotting = nanmean(cvs, 1);
    cvs_for_plotting(isnan(cvs_for_plotting)) = 0;

end
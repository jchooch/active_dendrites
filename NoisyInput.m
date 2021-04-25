%% NoisyInput.m

function [currents, muvec, movave] = NoisyInput(mus, tvec)
    
    dt = 0.1;
    tau = 3;
    sigma = 0.3;
    currents = [];
    muvec = [];
    for mu = mus
        current = [];
        current(1) = mu;
        for i = 2:(length(tvec) / length(mus))
            current(i) = current(i-1) + ((mu - current(i-1)) / tau)*dt + ...
                sigma*normrnd(0,1)*sqrt((2*dt) / tau);
        end
        currents = horzcat(currents, current);
        muv = ones(size(current)) .* mu;
        muvec = horzcat(muvec,muv);
    end
    for i = length(currents):length(tvec)
        currents(i) = mu(end);
    end
    movave = movmean(currents, 1000);
    
end
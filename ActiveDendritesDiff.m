%%writefile ActiveDendritesDiff.m

function ActiveDendritesDiff(t, y, C) %tspan, yinit, constants

dydt = [y1;
        y2;
        y3;
        y4]

end


% I(t+dt) = I(t) + ((mu - I(t))/tau)*dt + (sigma * G_t * sqrt(2*d*t/tau))
% dI/dt = ((mu - I)/tau)*dt + (sigma * G_t * sqrt(2*d*t/tau))
% dV_S/dt = ((1/R_S) * (V_rest_S - V_s) + (1/R_T) * (V_D - V_S) + I_AHP +
% I_app_S) / C_S
% dV_D/dt = ((1/R_D) * (V_rest_D - V_D) + (1/R_T) * (V_S - V_D) + I_Ca +
% I_app_D) / C_D
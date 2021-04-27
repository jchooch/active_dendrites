%function hinf

function output = hinf(v)
    output = -0.5*v -10;
    if output > 1
        output = 1;
    end
    if output < 0
        output = 0;
    end
    %output(v<=-22) = 1;
    %output(v>=-20) = 0;
end
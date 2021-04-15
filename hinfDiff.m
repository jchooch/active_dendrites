%function hinf

function output = hinf(v)
    output = -0.5*v -10;
    output(v<=-22) = 1;
    output(v>=-20) = 0;
end
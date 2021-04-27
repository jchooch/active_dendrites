%function minf

function output = minf(v)
    output = 0.5*v +5;
    if output > 1
        output = 1;
    end
    if output < 0
        output = 0;
    end
    %output(v<=-10) = 0;
    %output(v>=-8) = 1;
end
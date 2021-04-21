%function minf

function output = minf(v)
    output = 0.5*v +5;
    output(v<=-10) = 0;
    output(v>=-8) = 1;
end
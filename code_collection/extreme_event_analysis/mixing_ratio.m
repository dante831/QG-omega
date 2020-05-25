function r = mixing_ratio(vapor_pressure, pressure)

    rdgas = 287.04;
    rvgas = 461.50;

    r = rdgas * vapor_pressure ./ rvgas ./ (pressure - vapor_pressure);

end

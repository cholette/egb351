function pax = wind_rose(directions)

polarhistogram(pi/180*directions,'BinWidth',pi/180*10, ...
    'Normalization','probability');
pax = gca;
pax.RTick = linspace(pax.RTick(1),pax.RTick(end),5);
pax.ThetaZeroLocation = 'top'; % wind roses have N at top
pax.ThetaDir = "clockwise";    % wind roses have + clockwise

idx = find(pax.ThetaTick==0);
pax.ThetaTickLabel{idx} = pax.ThetaTickLabel{idx}+" (N)";

idx = find(pax.ThetaTick==90);
pax.ThetaTickLabel{idx} = pax.ThetaTickLabel{idx}+" (E)";

idx = find(pax.ThetaTick==180);
pax.ThetaTickLabel{idx} = pax.ThetaTickLabel{idx}+" (S)";

idx = find(pax.ThetaTick==270);
pax.ThetaTickLabel{idx} = pax.ThetaTickLabel{idx}+" (W)";
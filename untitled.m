

T = eccGCL(ii) - xi;
Displacement = ((ai/6 .* T + bi/2) .* T + ci) .* T + di ;

eccIS(ii) = eccGCL(ii) - Displacement(eccRange);
disp_dist(ii) = Displacement(eccRange);
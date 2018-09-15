function DisplacemantTP_Zeiss


%%
tp  = readtable('10-2testpoint.csv');

%% Drasdo model
% Drasdo N. "The length of Henle fibers in the human retina and a model of
% ganglion receptive field density in the visual field."
% Vison Research 2007.
%
% Lateral displacement at a location in the ganglion cell layer (GCL)
% is calculated using the coefficients ai, bi, ci, and di
% and a temporary variable T as follows:
% for an eccentricity in the ganglion cell layer (eccGCL) falling
% between xi and xi + 1,


%% Nasal [mm]
xi  = [0, 0.6243, 2.6231];
xi1 = [0.6243, 2.6231, 3.9632];

ai = [-4.3774, 1.2022, 0];
bi = [1.1856, -1.5470, 0];
ci = [0.6898, 0.5770, -0.1098];
di = [0, 0.4841, 0.147];

%% displacement from gcc leyer
eccGCL = 0.1 : 0.01 : 5 ;%: 30;
for ii  = 1 : length(eccGCL)
    
    % piecewise function depending on eccentricity
    if xi(1) <= eccGCL(ii) && xi(2) >= eccGCL(ii)
        eccRange = 1;
    elseif xi(2) <= eccGCL(ii) && xi(3) >= eccGCL(ii)
        eccRange = 2;
    elseif xi(3) <= eccGCL(ii)
        eccRange = 3;
    end
    
    T = eccGCL(ii) - xi;
    Displacement = ((ai/6 .* T + bi/2) .* T + ci) .* T + di ;
    
    eccIS(ii) = eccGCL(ii) - Displacement(eccRange);
    disp_dist(ii) = Displacement(eccRange);
    
end

% disp_dist > 0

figure; hold on;
plot( eccIS(disp_dist > 0),  disp_dist(disp_dist > 0), '-')
%% Temporal
xi  = [0, 1.2337, 2.5360]; % eccentricity range
xi1 = [1.2337, 2.5360, 5];

ai = [-0.103, 1.3537, 0];
bi = [-0.765, -0.8921, 0];
ci = [0.9336, -0.0885, -0.0689];
di = [0, 0.5374, 0.1639];

% displacement from gcc leyer
eccGCL = 0.1 : 0.1 : 5 ;%: 30;
for ii  = 1 : length(eccGCL)
    
    % piecewise function depending on eccentricity
    if xi(1) <= eccGCL(ii) && xi(2) >= eccGCL(ii)
        eccRange = 1;
    elseif xi(2) <= eccGCL(ii) && xi(3) >= eccGCL(ii)
        eccRange = 2;
    elseif xi(3) <= eccGCL(ii)
        eccRange = 3;
    end
    
    T = eccGCL(ii) - xi;
    Displacement = ((ai/6 .* T + bi/2) .* T + ci) .* T + di ;
    
    eccIS(ii) = eccGCL(ii) - Displacement(eccRange);
    disp_dist(ii) = Displacement(eccRange);
    
end

% add temporal
plot( eccIS(disp_dist > 0),  disp_dist(disp_dist > 0), '-')
title('lateral displacement in mm')

%% add Sup and Inf displacement
plot( eccIS,  disp_dist * 0.76, '-')

legend({'Nasal', 'Temporal', 'Sup/Inf'})

ax = gca;
ax.XAxisLocation ='origin';
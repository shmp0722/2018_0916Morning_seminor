function DiplacementTP

%% 10-2 test point locations. convensional and displaced.
% tp_disp = readtable('10-2testpoint_displacement.xlsx');
tp  = readtable('10-2testpoint.csv');

% %% eccentricity
% tp_disp.ecc(abs(tp_disp.x)>abs(tp_disp.y)) = abs(tp_disp.x(abs(tp_disp.x)>abs(tp_disp.y)));
% tp_disp.ecc(abs(tp_disp.x)<abs(tp_disp.y)) = abs(tp_disp.y(abs(tp_disp.x)<abs(tp_disp.y)));
%

%% conventinal test point

figure; hold on;


% add circle
R = [1, 3, 5, 7, 9]; % radious
C = jet(length(R));  % color for lines

cx = 0; cy = 0; % center

t = linspace(0,2*pi,100);

for i = 1: length(R)
    r = R(i);           % ??
    plot(r*sin(t)+cx,r*cos(t)+cy,'Color',C(i,:), 'LineWidth',2.5)
end

legend(num2str(R(1)),num2str(R(2)),num2str(R(3)),num2str(R(4)),num2str(R(5)))


plot(tp.x, tp.y,'sk','MarkerSize',10)%, 'MarkerFaceColor','k');
axis equal

axis square
title 'Conventinal test point'
set(gca, 'FontSize',18)
%%
saveas(gca, fullfile(pwd,'/Figure','ConventionalTestPoint.png'))

%% Sj?strand J. Graefe?s Arch Clin Exp Ophthalmol 1999
% x = Cone ecc [mm]
% X = x/3.6 [degree]

%
% M_angle = atan2(tp.y(i),tp.x(i))*180/pi; % sita = atan2(Y,X)
%
% M_angle = atan2(tp.y(i),tp.x(i)); % sita = atan2(Y,X)

% disp_y = 0.37*exp(-((x-0.67)/1.12)^2);

%% Sjostrand formula

disp_mm = 1.29*(tp.ecc+0.046).^0.67; %in [mm]

% disp_deg = disp_mm./3.6; % convert mm in deg
disp_deg = disp_mm./2.86; % Cirrus assumption convert mm in deg


tp.disp_mm  = disp_mm; % distance displacement
tp.disp_deg = disp_deg; % convert to deg

tp.Theta =  atan2(tp.y,tp.x); % angle of each test point

tp.disp_x = (tp.ecc+disp_deg) .* cos(tp.Theta); %
tp.disp_y = (tp.ecc+disp_deg) .* sin(tp.Theta); %

%% make figure to show this Sjostrand's model
figure; hold on;

% add circle
R = [1, 3, 5, 7, 9];

C = jet(length(R));

cx = 0; cy = 0; % ??

t = linspace(0,2*pi,100);

for i = 1: length(R)
    r = R(i);           % ??
    plot(r*sin(t)+cx,r*cos(t)+cy,'Color',C(i,:), 'LineWidth',2.5)
end

legend(num2str(R(1)),num2str(R(2)),num2str(R(3)),num2str(R(4)),num2str(R(5)))

plot(tp.disp_x, tp.disp_y,'or','MarkerSize',10)%, 'MarkerFaceColor','k');
plot(tp.x, tp.y,'sk','MarkerSize',8)%, 'MarkerFaceColor','k');

axis equal
title 'Sjostrand model'
set(gca,'FontSize',18)
%% make figure to show this Sjostrand's model
figure; hold on;

% add circle
R = [1, 3, 5, 7, 9];

C = jet(length(R));

cx = 0; cy = 0; % ??

t = linspace(0,2*pi,100);

for i = 1: length(R)
    r = R(i);           % ??
    plot(r*sin(t)+cx,r*cos(t)+cy,'Color',C(i,:), 'LineWidth',2.5)
end

legend(num2str(R(1)),num2str(R(2)),num2str(R(3)),num2str(R(4)),num2str(R(5)))

plot(tp.disp_x, tp.disp_y,'sk','MarkerSize',10)%, 'MarkerFaceColor','k');
% plot(tp.x, tp.y,'sk','MarkerSize',8)%, 'MarkerFaceColor','k');

axis equal
title 'Sjostrand model'
set(gca,'FontSize',18)

%%
saveas(gca, fullfile(pwd,'Figure/SjostrandModel.png'))

%% Sjostrand J. Graefe?s Arch Clin Exp Ophthalmol 1999
% x = Cone ecc [mm]
% X = x/3.6 [degree]


M_angle = atan2(tp.y(i),tp.x(i));% *180/pi; % sita = atan2(Y,X)

ecc_deg = sqrt(tp.y(i)^2+tp.x(i)^2);
ecc_mm  = 3.6*ecc_deg;

displ_mm = 0.37*exp(-((ecc_mm-0.67)/1.12)^2);
displ_deg = displ_mm/3.6;


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

% add Sup and Inf displacement
plot( eccIS,  disp_dist * 0.76, '-')

legend({'Nasal', 'Temporal', 'Sup/Inf'})

ax = gca;
ax.XAxisLocation ='origin';


%% Trupin described Drasdo's formula like this
% 1. Dc(r, Theta) ; density of RGCs abd Df(r,Theta) density of RGC
% receptive fields at some retinal location(r, Theta) in polar coordinates.
% Btoh functions give density in [cells/mm^2]

% 2. Assume we wish to estimate the displacement for receptive field
% location (R,h), and that alpha(degree) and K are chosen so that the sector subtending
% alpha degrees centered on h is tiled into K subsectors (as in the left panel
% of Figure 1).

% 3. Compute the number of receptive fields in the sector out to
% radius R mm:
%
% Nf = symsum((alpha*pi/360)*(2*iR^2/K^2)*Df(iR/K, theta), k , 1, K)
%
% 4.Compute the number of RGCs in the sector out to radius kR/K mm,
% for some k >= 1:
%
% Nc = symsum((alpha*pi/360)*(2*iR^2/K^2)*Dc(iR/K, theta), k , 1, K)
%
% This is the RGC curve in the right of Figure 1.
%
% 5. Find k so that jNc(k)  Nfj is minimized.
%
% 6. Report displacement as kR/K  R mm.
%
% In this paper, we used the equations of Drasdo et al.6 to
% compute Df, at angular eccentricity (E) as follows:
%
% Rve = Rv0(1 + r/E2v(theta));
%
% Roe = Ro0(1 + r/E2o(theta));
%
% k(r) = 1+(1.004 - 0.007209*r + 0.001694*r^2 - 0:00003765^r^3)^-2;
%
% Df(r/theta) = (1.12 + 0.0273*r)*k/1.155/(Rve^2-Roe^2);

%%
% tp  = readtable('10-2testpoint.csv');
tp.Theta =  atan2(tp.y,tp.x); % angle of each test point

tp.turpin_disp_x = (tp.ecc+tp.turpin_disp) .* cos(tp.Theta); %
tp.turpin_disp_y = (tp.ecc+tp.turpin_disp) .* sin(tp.Theta); %

%% save progression
% writetable(tp,'10-2testpoint.csv')
%% figure
figure; hold on;

% add circle
R = [1, 3, 5, 7, 9];
% R = [3.4, 5.6, 6.8, 8.3, 9.7];

C = jet(length(R));

cx = 0; cy = 0; % ??

t = linspace(0,2*pi,100);

for i = 1: length(R)
    r = R(i);           % ??
    plot(r*sin(t)+cx,r*cos(t)+cy,'Color',C(i,:), 'LineWidth',2.5)
end

legend(num2str(R(1)),num2str(R(2)),num2str(R(3)),num2str(R(4)),num2str(R(5)))

plot(tp.turpin_disp_x , tp.turpin_disp_y,'sk','MarkerSize',10)%, 'MarkerFaceColor','k');
% plot(tp.x, tp.y,'sk','MarkerSize',8)%, 'MarkerFaceColor','k');

axis equal
title 'Drasdo-Turpin model'
set(gca,'FontSize',18)

%% save the figure
saveas(gca, fullfile(pwd,'/Figure/Drasdo-TurpinModel.png'))

%% 10-2 turpin figure
figure; 
subplot(1,2,1)
hold on;
plot(tp.x , tp.y,'sk','MarkerSize',10)%, 'MarkerFaceColor','k');

axis equal
title 'conventional'
set(gca,'FontSize',18)

subplot(1,2,2)
hold on;
plot(tp.turpin_disp_x , tp.turpin_disp_y,'sk','MarkerSize',10, 'MarkerFaceColor','r');
% plot(tp.x, tp.y,'sk','MarkerSize',8)%, 'MarkerFaceColor','k');

axis equal
title 'Drasdo-Turpin model'
set(gca,'FontSize',18)
get(gca,'XTickLabel')

ax = gca;
% ax.XAxisLocation = 'origin';

%% 24-2c
tp_new  = readtable('24-2cXYcoordinates.xlsx','Sheet',2);

%% 24-2c turpin figure
figure; 
subplot(1,2,1)
hold on;
plot(tp.x , tp.y,'sk','MarkerSize',10)%, 'MarkerFaceColor','k');

plot(tp_new.new_x , tp_new.new_y,'sr','MarkerSize',10, 'MarkerFaceColor','r');

axis equal
title 'conventional'
set(gca,'FontSize',18)

subplot(1,2,2)
hold on;
plot(tp.turpin_disp_x , tp.turpin_disp_y,'sk','MarkerSize',10)%, 'MarkerFaceColor','r');
plot(tp_new.new_turpin_disp_x(1:10) , tp_new.new_turpin_disp_y(1:10), ...
    'sr','MarkerSize',10, 'MarkerFaceColor','r');


% plot(tp.x, tp.y,'sk','MarkerSize',8)%, 'MarkerFaceColor','k');

axis equal
title 'Drasdo-Turpin model'
set(gca,'FontSize',18)
get(gca,'XTickLabel')

ax = gca;
% ax.XAxisLocation = 'origin';

%%
for ii = 1 :  sum(~isnan( tp_new.new_x))
    %finding same testpoint
    a =  tp.x == tp_new.new_x(ii); 
    b =  tp.y == tp_new.new_y(ii);
    
    tp_new.new_turpin_disp_x(ii) = tp.turpin_disp_x(a + b == 2) ;
    tp_new.new_turpin_disp_y(ii) = tp.turpin_disp_y(a + b == 2) ;
    
    tp_new.new_turpin_disp(ii) = tp.turpin_disp(a + b == 2) ;
end

%%



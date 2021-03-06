function TestPoint
%% load applications

% addpath(genpath('spe')); % subplot package
% addpath(genpath('mtit')); % subplot package

%% dipict 10-2 test point and 10, 20 and 30 eccentricities.

tp_10 = readtable('10-2testpoint.csv');
tp_30 = readtable('30-2testpoint.xlsx');
tp_24 = readtable('24-2cXYcoordinates.xlsx', 'Sheet',2);

%% 10-2 only
figure; hold on;
% for ii = 1: length(tp_10.x)
    plot(tp_10.x, tp_10.y, 'sk')%, 'MarkerFaceColor', c(m(ii),:),...
    %             'MarkerEdgeColor', c(m(ii),:),'MarkerSize',12)
% end

% t = linspace(-pi, pi, 100);
% plot(10*sin(t),10*cos(t))
% plot(20*sin(t),20*cos(t))
% plot(30*sin(t),30*cos(t))

axis equal
axis on
set(gca, 'XLim',[-11 11],'YLim',[-11 11])
%     grid on
%     Grid = -12:3:12;
%     set(gca,'XTick',Grid)
%     set(gca,'YTick',Grid)
saveas(gca, 'figure/10-2.png')
%% 10-2 new only
figure; hold on;
% for ii = 1: length(tp_10.x)
    plot(tp_24.conv_24_x, tp_24.conv_24_y, 'sk','MarkerSize',10)%, 'MarkerFaceColor', c(m(ii),:),...
    %             'MarkerEdgeColor', c(m(ii),:),'MarkerSize',12)
% end
 plot(tp_24.new_x, tp_24.new_y, 'sk','MarkerSize',10,...
      'MarkerFaceColor','r')

% t = linspace(-pi, pi, 100);
% plot(10*sin(t),10*cos(t))
% plot(20*sin(t),20*cos(t))
% plot(30*sin(t),30*cos(t))

axis equal
axis on
set(gca, 'XLim', [-35,35])
set(gca, 'YLim', [-35,35])
set(gca, 'XTick', -30:10:30)
set(gca, 'YTick', -30:10:30)

ax = gca;
ax.XAxisLocation = 'origin'
ax.YAxisLocation = 'origin'


%% 10-2 & 30-2 test point

figure; hold on;
% dipic 30-2
    plot(tp_30.x, tp_30.y, 'sk','MarkerSize',12)%, 'MarkerFaceColor', c(m(ii),:),...
    %             'MarkerEdgeColor', c(m(ii),:),'MarkerSize',12)

% dipic 10-2
    plot(tp_10.x, tp_10.y, 'sr')%, 'MarkerFaceColor', c(m(ii),:),...
    %             'MarkerEdgeColor', c(m(ii),:),'MarkerSize',12)

% 10, 20, 30 degree circles
% t = linspace(-pi, pi, 100);
% plot(10*sin(t),10*cos(t))
% plot(20*sin(t),20*cos(t))
% plot(30*sin(t),30*cos(t))

% add 8*8 spectralis grid
% Grid = -12:3:12;
% % figure; hold on;
% for ii = 1 : length(Grid)
%     plot([Grid(1), Grid(end)],[Grid(ii), Grid(ii)],'-k')
%     plot([Grid(ii), Grid(ii)],[Grid(1), Grid(end)],'-k')
% end

axis equal
% axis setting
set(gca, 'XLim', [-35,35])
set(gca, 'YLim', [-35,35])
set(gca, 'XTick', -30:10:30)
set(gca, 'YTick', -30:10:30)
%% displaced 10-2 & 30-2 test point

% Sjostrand formula
tp_30.ecc = sqrt( tp_30.x.^2 + tp_30.y.^2);

disp_mm = 1.29*(tp_30.ecc+0.046).^0.67; %in [mm]

% disp_deg = disp_mm./3.6; % convert mm in deg 
disp_deg = disp_mm./  3.4965; % convert mm in deg ; Cirrus assumption


tp_30.disp_mm  = disp_mm; % distance displacement
tp_30.disp_deg = disp_deg; % convert to deg

tp_30.Theta =  atan2(tp_30.y,tp_30.x); % angle of each test point

tp_30.disp_x = (tp_30.ecc+disp_deg) .* cos(tp_30.Theta); % 
tp_30.disp_y = (tp_30.ecc+disp_deg) .* sin(tp_30.Theta); % 

%% displaced 10-2 & 30-2 test point
figure; 
subplot(1,2,2); hold on;
% dipic 30-2
    plot(tp_30.disp_x, tp_30.disp_y, 'sk','MarkerSize',12)%, 'MarkerFaceColor', c(m(ii),:),...
    %             'MarkerEdgeColor', c(m(ii),:),'MarkerSize',12)

% dipic 10-2
    plot(tp_10.turpin_disp_x, tp_10.turpin_disp_y, 'sr')%, 'MarkerFaceColor', c(m(ii),:),...
    %             'MarkerEdgeColor', c(m(ii),:),'MarkerSize',12)

% % add 8*8 spectralis grid
% Grid = -12:3:12;
% % figure; hold on;
% for ii = 1 : length(Grid)
%     plot([Grid(1), Grid(end)],[Grid(ii), Grid(ii)],'-k')
%     plot([Grid(ii), Grid(ii)],[Grid(1), Grid(end)],'-k')
% end

axis equal
% axis setting
set(gca, 'XLim', [-35,35])
set(gca, 'YLim', [-35,35])
set(gca, 'XTick', -30:10:30)
set(gca, 'YTick', -30:10:30)


title('Displacement')


% 10-2 & 30-2 test point

subplot(1,2,1); hold on;
% dipic 30-2
    plot(tp_30.x, tp_30.y, 'sk','MarkerSize',12)%, 'MarkerFaceColor', c(m(ii),:),...
    %             'MarkerEdgeColor', c(m(ii),:),'MarkerSize',12)

% dipic 10-2
    plot(tp_10.x, tp_10.y, 'sr')%, 'MarkerFaceColor', c(m(ii),:),...
    %             'MarkerEdgeColor', c(m(ii),:),'MarkerSize',12)

% 10, 20, 30 degree circles
% t = linspace(-pi, pi, 100);
% plot(10*sin(t),10*cos(t))
% plot(20*sin(t),20*cos(t))
% plot(30*sin(t),30*cos(t))

% % add 8*8 spectralis grid
% Grid = -12:3:12;
% % figure; hold on;
% for ii = 1 : length(Grid)
%     plot([Grid(1), Grid(end)],[Grid(ii), Grid(ii)],'-k')
%     plot([Grid(ii), Grid(ii)],[Grid(1), Grid(end)],'-k')
% end

axis equal
% axis setting
set(gca, 'XLim', [-35,35])
set(gca, 'YLim', [-35,35])
set(gca, 'XTick', -30:10:30)
set(gca, 'YTick', -30:10:30)

title('Conventional')
% mtit('HAF test point')
%%
saveas(gca, 'GridTestPoints.png')
function HFA30_2vs24_2

%% load
tp_new = readtable('24-2cXYcoordinates.xlsx','Sheet',2);
% tp_disp = readtable('10-2testpoint_displacement.xlsx');
tp  = readtable('10-2testpoint.csv');
%% 24-2
figure; hold on;

% add circle
% R = [1, 3, 5, 7, 9]; % radious
% C = jet(length(R));  % color for lines
% 
% cx = 0; cy = 0; % center
% 
% t = linspace(0,2*pi,100);
% 
% for i = 1: length(R) 
%     r = R(i);           % ??
%     plot(r*sin(t)+cx,r*cos(t)+cy,'Color',C(i,:), 'LineWidth',2.5)
% end
% legend(num2str(R(1)),num2str(R(2)),num2str(R(3)),num2str(R(4)),num2str(R(5)))

plot(tp_new.conv_24_x, tp_new.conv_24_y,'sk','MarkerSize',10)%, 'MarkerFaceColor','k');

title 'Conventional 24-2'
set(gca, 'FontSize',18)

axis equal
axis square
set(gca,'XLim',[-30 30], 'YLim',[-30 30])

% plot(tp_new.new_x, tp_new.new_y,'sr','MarkerSize',10, 'MarkerFaceColor','r'...
%     ,'MarkerEdgeColor', 'k')%, 'MarkerFaceColor','k');

set(gca, 'XTick',[-30:10:30],'FontSize',10)

ax = gca ;
ax.XAxisLocation   
ax.XAxisLocation   = 'origin'; 
ax.YAxisLocation   = 'origin'; 

%%
saveas(gca, '24-2.pdf')
saveas(gca, '24-2.png')

%% 24-2c
figure; hold on;

% add circle
% R = [1, 3, 5, 7, 9]; % radious
% C = jet(length(R));  % color for lines
% 
% cx = 0; cy = 0; % center
% 
% t = linspace(0,2*pi,100);
% 
% for i = 1: length(R) 
%     r = R(i);           % ??
%     plot(r*sin(t)+cx,r*cos(t)+cy,'Color',C(i,:), 'LineWidth',2.5)
% end
% legend(num2str(R(1)),num2str(R(2)),num2str(R(3)),num2str(R(4)),num2str(R(5)))

plot(tp_new.conv_24_x, tp_new.conv_24_y,'sk','MarkerSize',10)%, 'MarkerFaceColor','k');

title 'Conventional 24-2'
set(gca, 'FontSize',18)

axis equal
axis square
set(gca,'XLim',[-30 30], 'YLim',[-30 30])

plot(tp_new.new_x, tp_new.new_y,'sr','MarkerSize',10, 'MarkerFaceColor','r'...
    ,'MarkerEdgeColor', 'k')%, 'MarkerFaceColor','k');

set(gca, 'XTick',[-30:10:30],'FontSize',10)

ax = gca ;
ax.XAxisLocation   
ax.XAxisLocation   = 'origin'; 
ax.YAxisLocation   = 'origin'; 

%%
saveas(gca, '24-2c.pdf')
saveas(gca, '24-2c.png')

%% 24-2c
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

plot(tp_new.conv_24_x, tp_new.conv_24_y,'sk','MarkerSize',10)%, 'MarkerFaceColor','k');

title 'Conventional 24-2'
set(gca, 'FontSize',18)

axis equal
axis square
set(gca,'XLim',[-30 30], 'YLim',[-30 30])

plot(tp_new.new_x, tp_new.new_y,'sr','MarkerSize',10, 'MarkerFaceColor','r'...
    ,'MarkerEdgeColor', 'k')%, 'MarkerFaceColor','k');

set(gca, 'XTick',[-30:10:30],'FontSize',10)

ax = gca ;
ax.XAxisLocation   
ax.XAxisLocation   = 'origin'; 
ax.YAxisLocation   = 'origin'; 

%%
saveas(gca, '24-2cWithEcc.pdf')
saveas(gca, '24-2cWithEcc.png')

%% Sjostrand J. Graefe?s Arch Clin Exp Ophthalmol 1999 
% x = Cone ecc [mm] 
% X = x/3.6 [degree]


M_angle = atan2(tp_24.new_y,tp_24.new_x);% *180/pi; % sita = atan2(Y,X)

ecc_deg = sqrt(tp_24.new_y .^2 + tp_24.new_x .^2);
ecc_mm  = 3.6*ecc_deg;

displ_mm = 0.37*exp(-((ecc_mm-0.67)/1.12).^2);
displ_deg = displ_mm/3.6;





%%
tp_new.Theta_all =  atan2(tp_new.x24_2c_x,tp_new.x24_2c_y); % angle of each test point
tp_new.Theta_new =  atan2(tp_new.new_x,tp_new.new_y); % angle of each test point


tp_new.ecc_all   =  sqrt( tp_new.x24_2c_x .^2 + tp_new.x24_2c_y .^2);
tp_new.ecc_new   =  sqrt( tp_new.new_x .^2 + tp_new.new_y .^2);

%% eccentricity 
tp_new.ecc(abs(tp_new.new_x)>abs(tp_new.new_y)) = abs(tp_new.new_x(abs(tp_new.new_x)>abs(tp_new.new_y)));
tp_new.ecc(abs(tp_new.new_x)<abs(tp_new.new_y)) = abs(tp_new.new_y(abs(tp_new.new_x)<abs(tp_new.new_y)));
 
tp_new.turpin_disp_x = (tp_new.ecc+tp_new.turpin_disp) .* cos(tp_new.Theta); % 
tp_new.turpin_disp_y = (tp_new.ecc+tp_new.turpin_disp) .* sin(tp_new.Theta); % 



%% figure

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


plot(tp_new.new_x, tp_new.new_y,'sk','MarkerSize',10)%, 'MarkerFaceColor','k');
axis equal

axis square
title 'Newly added test point to 24-2C'
set(gca, 'FontSize',18)

%%
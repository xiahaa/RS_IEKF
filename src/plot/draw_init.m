clc;close all;
load('draw_res.mat');

% set the initial values of parameters
para.ts = 0.025; % readout time of rolling shutter camera, 0.025 
para.td = 1.175;%1.175; % timestamp delay	
para.wd = [0 0 0]; % gyroscope bias
para.f = 700; % camera focal length
para.fx = 671; % camera focal length
para.fy = 671; % camera focal length
para.w = 720; % frame width (this is a fixed value)
para.h = 480; % frame height (this is a fixed value)
para.cx = (para.w-1)/2; % principal point
para.cy = (para.h-1)/2; % principal point
para.rcam = [pi/sqrt(2); -pi/sqrt(2); 0]; % relative orientation between camera and gyro
para.sigma = 1e-6; % the noise variance of gyro readings
para.pn = 1; % noise variance of feature detection

sigma_draw = [20,20,20,0.1,0.01];

nomial.f = para.f;
nomial.cx = para.cx;
nomial.cy = para.cy;
nomial.ts = para.ts;
nomial.td = para.td;
nomials = [nomial.f,nomial.cx,nomial.cy,nomial.ts,nomial.td];

valid_res1 = res1;
valid_res2 = res2;
valid_res3 = res3;

font_size = 15;
num_bins = 8;
subplot_margin = 0.1;
subplot_spacing = 0.08;

cmap = lines(3);

x0 = 0;
y0 = 0;
width = 8;
height = 10;
fig = figure('Units','inches',...
'Position',[x0 y0 width height]);

ind = [1,2,3,4,5];

tlname = {'f: pixel', 'c_u: pixel', 'c_v: pixel', 't_r: (s)', 't_d: (s)'};

for i = 1:length(ind)
    ii = ind(i);
    subaxis(3,2,i, 'Margin', subplot_margin, 'Spacing', subplot_spacing);
    histogram(valid_res1(:,ii),linspace(0,nomials(ii)*2,num_bins), 'Normalization', 'Probability', ...
    'EdgeColor', cmap(1,:), 'LineWidth', 0.5, 'FaceColor', cmap(1,:), 'FaceAlpha', 0.2)
    hold on;
    histogram(valid_res2(:,ii),linspace(0,nomials(ii)*2,num_bins), 'Normalization', 'Probability', ...
            'EdgeColor', cmap(2,:), 'LineWidth', 0.5, 'FaceColor', cmap(2,:), 'FaceAlpha', 0.4)
    histogram(valid_res3(:,ii),linspace(0,nomials(ii)*2,num_bins), 'Normalization', 'Probability', ...
            'EdgeColor', cmap(3,:), 'LineWidth', 0.5, 'FaceColor', cmap(3,:), 'FaceAlpha', 0.6)
    plot([nomials(ii),nomials(ii)],get(gca,'YLim'),'r--','LineWidth',2);
    hold off;
    grid minor
    set(gca,'FontSize', font_size);
    %xticklabels({round(linspace(0,nomial.f*2,num_bins))})
    %xlabel('f: pixel','FontSize', font_size, 'Interpreter', 'latex')
    title(tlname{ii},'FontSize', font_size)
    
    if i == 1
        lgnd = legend({'[13]','[14]','Proposed','Init'}, 'Location', 'NorthEast');
        set(lgnd,'FontSize', font_size-2);
    end
end

subaxis(3,2,6, 'Margin', subplot_margin, 'Spacing', subplot_spacing);
bar(1,sum(failures1==0)/length(failures1),'EdgeColor', cmap(1,:), 'LineWidth', ...
    0.5, 'FaceColor', cmap(1,:), 'FaceAlpha', 0.2);
hold on;
bar(2,sum(failures2==0)/length(failures2),'EdgeColor', cmap(2,:), 'LineWidth', ...
    0.5, 'FaceColor', cmap(2,:), 'FaceAlpha', 0.4);
bar(3,sum(failures3==0)/length(failures3),'EdgeColor', cmap(3,:), 'LineWidth', ...
    0.5, 'FaceColor', cmap(3,:), 'FaceAlpha', 0.6);
xticks(1:3);
xticklabels({'[13]','[14]','Proposed'});
xtickangle(20);
title('Success rate');
grid minor;
hold off;
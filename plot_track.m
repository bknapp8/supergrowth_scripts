% Plot a single track
close all;
k = randi(length(TRACKS));
% k = 419;
res = 0.08;

x = TRACKS(k).x;
y = TRACKS(k).y;

figure;
plot(x*res,y*res,'b-');
axis equal;
box off
xm = mean(x)*res;
ym = mean(y)*res;
xlim([xm-1 xm+1]);
ylim([ym-1 ym+1]);


%%
dt = TRACKS(k).frame_interval;

% Use function from Holt lab for MSD
[MSD,MSD_err] = calculate_MSD_err(x,y,0,dt,res);
% Make time vector
time = (0:1:(length(MSD)-1))'*dt;

figure;
plot(time, MSD, 'r');
xlim([0 1000]);
ylim([0 1]);
set(gca, 'fontsize', 20);
xlabel('Time (ms)');
ylabel('Mean squared displacement ()');


% Errorbar plot from MSD calculation
figure;
errorbar(time,MSD,MSD_err, 'ro-');
set(gca, 'fontsize', 20);
xlabel('dt (ms)');
ylabel('Mean squared displacement (\mum^2)');

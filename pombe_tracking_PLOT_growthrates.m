close all;

t_int = 5; % Frame interval in minutes
res = 0.11; % Pixel resolution in microns/pixel

tree = TOTALS;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Vectors for storing all data points
bigT = [];
bigR = [];
bigL = [];
figure;

% For counting cell trajectories of proper length
Ncell = 0;

for k=1:length(tree)
    clearvars L T;
    endcut = 8; % % Use 0 if no cutting at end
    begcut =5; % Use 1 if no cutting at beginning
    if length(tree(k).traj)>20
        
        smf = 5; % Smoothing factor
        for j=1:length(tree(k).traj)
            L(j) = tree(k).traj(j).area;
            T(j) = tree(k).traj(j).frame_num-1;
        end
        colr = rand(1,3);
        
        
        R = smooth(smooth(gradient(log(L)), smf),smf)/t_int*60;

        % Filter for noisy growth rates
        Rcut1 = 1;
        Rcut2 = -0.3;
        mean(R);
        noiseR = std(R)/mean(R);
        minR = min(R);
        maxL = 20;
        if (max(R)<Rcut1 & min(R)>Rcut2 & mean(R)>-0.2 & noiseR<2)
            Ncell = Ncell + 1;
            colr = rand(1,3);
            bigL =[bigL; smooth(L(begcut:end-endcut),smf)*res];
            normL = 1;
            bigT = [bigT; (T(begcut:end-endcut))'*t_int];
            bigR = [bigR; R(begcut:end-endcut)];
            
            %% Plot cell lengths
            subplot(2,1,1);
            plot(T(begcut:end-endcut)*t_int/60, (L(begcut:end-endcut))*res, 'color', colr, 'linewidth', 1);
            set(gca, 'fontsize', 20);
            xlabel('Time (h)')
            ylabel('Cell length (\mum)');
            
            % Label trajectories for correction
            if true
               text(max(T(end-endcut)*t_int)/60,(L(end-endcut)),num2str(k));
            end
            hold on;

            %% Plot growth rates
            subplot(2,1,2);   
            plot(T(begcut:end-endcut)*t_int/60, R(begcut:end-endcut),'color', colr, 'linewidth', 1);
            if true
               text(max(T(end-endcut)*t_int)/60,R(end-endcut),num2str(k));
            end
            set(gca, 'fontsize', 15);
            xlabel('Time (h)')
            ylabel('Size-Normalized Growth Rate (h^{-1})');
            hold on;
        end
            
    end
end

disp('Number of trajectories analyzed:');
disp(Ncell);


%% Plot binned rates with shaded error bars
gx = 0:t_int*5:max(bigT);
[b,n,s] = bindata(bigT,bigR, gx);

figure;
shadedErrorBar(gx,b,s, 'b', 0.1);
set(gca, 'fontsize', 20);
xlabel('Time (min)');
ylabel('Size-Normalized Growth Rate (h^{-1})');
set(gca, 'box', 'off');
ylim([-0.05 0.7]);
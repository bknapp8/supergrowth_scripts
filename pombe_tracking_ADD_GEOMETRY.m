% Script for adding cell geometry calculations to trajectories generated
% from pombe_tracking_MAKE_STRUCT and pombe_tracking_MAKE_TRAJECTORY


%% Loop through TOTALS structure to add geometries based on centerline calculation

stp = 10;

for k=1:length(TOTALS)
    for j=1:length(TOTALS(k).traj)
        x = TOTALS(k).traj(j).Xcont;
        y = TOTALS(k).traj(j).Ycont;
        
        [S,V,L,W] = CENTERLINE_geom_props(x,y,stp);
        
        TOTALS(k).traj(j).surface_area = S;
        TOTALS(k).traj(j).volume = V;
        TOTALS(k).traj(j).length = L;
        TOTALS(k).traj(j).width = W;
        
    end
end
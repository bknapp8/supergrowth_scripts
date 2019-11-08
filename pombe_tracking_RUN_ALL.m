% Script for analyzing CONTOURS data from S. pombe images. This script will
% generate a TOTALS structure with trajectories. This will also plot
% population-averaged growth rates over the timecourse of the experiment.

run pombe_tracking_MAKE_STRUCT;
run pombe_tracking_MAKE_TRAJECTORY;
run pombe_tracking_ADD_GEOMETRY;
run pombe_tracking_PLOT_growthrates;


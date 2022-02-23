%% SORTING
% Script Created by Connor Murphy
clc
clear all
close all

dumpfreq = 500; %IMPORTANT: determined by liggghts script
nsteps = 380000;  %IMPORTANT: determined by liggghts script
nt = nsteps/dumpfreq; %Number of timesteps
delta = 0.006562; %Pellet Diameter calculated in InertiaToRadius.m
timestep = 0.0001; %Seconds

dx = 4*delta;
dy = 2.5*delta;
dz = 2*delta;

name = 'dump.pellet_13_0D'; %Name of data file

maxP = 10; %Max amount of pellets that could fit in a bin (for preallocation) May want to increase for bigger bins
n1 = 12572; %IMPORTANT: Number of pellets in the first timestep. Find this value in dump file. Ex: ITEM: Number of Entries : Desired Value

[pInEachBin,COM,Vel, Quat, nx, ny, nz, FlowR]= binPelletAnalysisV2(name,nt,dx,dy,dz,maxP,n1); 
%Returns:
%pInEachBin: Pellet IDs telling which pellet is in which bin at a timestep
%pInEachBin: [xindex,yindex,zindex,timestep,pelletids]
%COM: Center of mass for each pellet
%COM: [Pellet,idxyz(4values),timestep]
%Vel: Velocity of center of mass for each pellet
%Vel: [Pellet,idxyz(4values),timestep]
%Quat: Quaternion values for each pellet
%Quat: [Pellet,idq0q1q2q3(5values),timestep]
%nx,ny,nz: Number of bins in each direction

%save('../Data_Plotting/PelletData2_0.mat','-v7.3') %Saving data to the Data_Plotting Directory. (version 7.3 for large files)
%NOTE: If you change the name of PelletData.mat, you have to add that name
%to the .gitignore file so it doesn't get pushed to GitHub

% FlowR13_0 = FlowR;
%save('PelletRate13_0.mat', 'FlowR13_0');






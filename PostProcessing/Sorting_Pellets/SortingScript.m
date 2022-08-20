%% SORTING
% Script Created by Connor Murphy
%***********************************************************
%This script reads data from a custom dump file produced by /scripts/in.fullbin
%(can use other scripts) that gives data in the following format: 
%id, quaternions, position, velocity: id q0 q1 q2 q3 xx xy xz vx vy vz
%It uses the pellet ids and positions for each timestep to create a 5D
%array(pInEAchbin) that stores pellets ids into spatial bins for each timestep
%Data is stored in the following way: [x, y, z, t, Pnum]
%It also takes data and creates 3D arrays for center of mass, quaternions,
%and velocity that store these components for each timestep. How they are
%stored is detailed in the function binPelletAnalysisV2 file.

%SETUP: 
%1)Make sure dumpfreq is set to the dump frequency of the dump file being
%used. This value can be found in the command that creates the dump file
%(currently every 500 timesteps)

%2)The number of timesteps needs to be set correctly, this is the total
%sum of all the values used for the run commands in the simulation script

%3)The timestep value needs to be correct. This value should be the same as
%the one used by the timestep command in the simulation script that was
%used.

%4)dx, dy, and dz are values that determine the size of the spatial bins
%pellets are sorted into. These can be refined or enlarged, but changing
%too drastically might cause issues that the user should be careful of. 

%5)Make sure that name is set to the correct dump file you want to use for
%sorting.

%6)Because pellets are added periodically to the simulation, this sorting
%script needs to know how many pellets there were for the first timestep.
%This is what n1 is used for. The value that needs to be used can be found
%at the beginning of the dump file. The dump file contains a line before
%each timestep that specifies how many pellets existed for that timestep. 

%***********************************************************

clc
clear all
close all

dumpfreq = 500; %IMPORTANT: determined by liggghts script
nsteps = 380000;  %IMPORTANT: determined by liggghts script
nt = nsteps/dumpfreq; %Number of timesteps
delta = 0.0075; %Rough template pellet diameter
len = 0.025;
timestep = 0.0001; %Seconds

dx = 1*len;
dy = 2.6*delta;
dz = 17.44*delta;

name = 'dump.pellet_8_0DF'; %Name of data file

maxP = 60; %Max amount of pellets that could fit in a bin (for preallocation) May want to increase for bigger bins
n1 = 12572; %IMPORTANT: Number of pellets in the first timestep. Find this value in dump file. Ex: ITEM: Number of Entries : Desired Value

[pInEachBin,COM,Vel, Quat, nx, ny, nz, FlowR]= binPelletAnalysisV2(name,nt,dx,dy,dz,maxP,n1); 

%Returns:
%pInEachBin: Pellet IDs telling which pellet is in which bin at each timestep
%pInEachBin: [xindex,yindex,zindex,timestep,pelletids]
%COM: Center of mass for each pellet
%COM: [Pellet,idxyz(4values),timestep]
%Vel: Velocity of center of mass for each pellet
%Vel: [Pellet,idxyz(4values),timestep]
%Quat: Quaternion values for each pellet
%Quat: [Pellet,idq0q1q2q3(5values),timestep]
%nx,ny,nz: Number of bins in each direction
%FlowR: The number of pellets pellets that have exited the hopper for
%each timestep. 

%save('../Data_Plotting/PelletData8_0F.mat','-v7.3') %Saving data to the Data_Plotting Directory. (version 7.3 for large files)
%NOTE: If you change the name of PelletData.mat, you have to add that name
%to the .gitignore file so it doesn't get pushed to GitHub

% FlowR5_0 = FlowR;
%save('PelletRate5_0.mat', 'FlowR5_0');








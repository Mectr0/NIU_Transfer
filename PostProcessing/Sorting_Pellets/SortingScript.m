%% SORTING
% Script Created by Connor Murphy
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
%pInEachBin: Pellet IDs telling which pellet is in which bin at a timestep
%pInEachBin: [xindex,yindex,zindex,timestep,pelletids]
%COM: Center of mass for each pellet
%COM: [Pellet,idxyz(4values),timestep]
%Vel: Velocity of center of mass for each pellet
%Vel: [Pellet,idxyz(4values),timestep]
%Quat: Quaternion values for each pellet
%Quat: [Pellet,idq0q1q2q3(5values),timestep]
%nx,ny,nz: Number of bins in each direction

%save('../Data_Plotting/PelletData8_0F.mat','-v7.3') %Saving data to the Data_Plotting Directory. (version 7.3 for large files)
%NOTE: If you change the name of PelletData.mat, you have to add that name
%to the .gitignore file so it doesn't get pushed to GitHub

% FlowR5_0 = FlowR;
%save('PelletRate5_0.mat', 'FlowR5_0');

%% Orientation confirmation

OrienVecz = zeros(1,3);
ConveyVec = [1 ; 0 ; 0];
z = [0; 0; 1];
for i = 1:nt
    for j=1:n1
    q0 = Quat(j,2,i);
    q1 = Quat(j,3,i);
    q2 = Quat(j,4,i);
    q3 = Quat(j,5,i);
    Ansz = [z(1)*(q0^2+q1^2-q2^2-q3^2)+2*z(2)*(q1*q2-q0*q3)+2*z(3)*((q0*q2)+(q1*q3));
            2*z(1)*(q0*q3+q1*q2) + z(2)*(q0^2-q1^2+q2^2-q3^2)+ 2*z(3)*((q2*q3)-(q0*q1));
            2*z(1)*(q1*q3-q0*q2) + 2*z(2)*(q0*q1 + q2*q3) + z(3)*(q0^2-q1^2-q2^2+q3^2)];
    OrienVecz(1,1)= Ansz(1,1);
    OrienVecz(1,2)= Ansz(2,1);
    OrienVecz(1,3)= Ansz(3,1);
    Gamma(1,j) = acosd((ConveyVec(1)*OrienVecz(1,1)+ ConveyVec(2)*OrienVecz(1,2) + ConveyVec(3)*OrienVecz(1,3)));
    end
    GammaMean(i) = mean(Gamma);
end

figure
plot ([1:80],GammaMean)
axis([0 80 -20 70])







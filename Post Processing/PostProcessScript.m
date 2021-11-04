%% SORTING
clear all
close all

dumpfreq = 1000; %determined by liggghts script
nsteps = 45000; %determined by liggghts script
nt = nsteps/dumpfreq; %Number of timesteps
delta = 0.006562; %Pellet Diameter 

dx = 4*delta;
dy = 2.5*delta;
dz = 2*delta;

name = 'dump.pellet_orienadd8'; %Name of data file

maxP = 15 ;%Max amount of pellets that could fit in a bin (for preallocation)
n1 = 10; %number of pellets in the first timestep

[pInEachBin,COM,Vel]= binPelletAnalysisV2(name,nt,dx,dy,dz,maxP,n1);

%% PLOTTING DATA

for i = 1:nt

Inspect = nonzeros(pInEachBin(:,1:50,:,i,:));
VelavgY(i) = mean(Vel(Inspect,3,i));

end
plot(1:nt,VelavgY,'-ro')
xlabel('Timestep')
ylabel('average Y-velocity (m/s)')

%(nx,ny,nz,nt,maxP)

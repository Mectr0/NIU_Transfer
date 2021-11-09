%% SORTING
clc
clear all
close all

dumpfreq = 1000; %determined by liggghts script
nsteps = 46000; %determined by liggghts script
nt = nsteps/dumpfreq; %Number of timesteps
delta = 0.006562; %Pellet Diameter 

dx = 4*delta;
dy = 2.5*delta;
dz = 2*delta;

name = 'dump.pellet_test'; %Name of data file

maxP = 8;%Max amount of pellets that could fit in a bin (for preallocation)
n1 = 9615; %number of pellets in the first timestep

[pInEachBin,COM,Vel,nx,ny,nz]= binPelletAnalysisV2(name,nt,dx,dy,dz,maxP,n1);

%% PLOTTING DATA

InspectXlo = 1;
InspectXhi = 19;
InspectYlo = 10;
InspectYhi = 61;
InspectZlo = 1;
InspectZhi = 10;

for i = 1:nt

Inspect = nonzeros(pInEachBin(:,10:60,:,i,:));
VelavgY(i) = mean(Vel(Inspect,3,i));
Inspect2 = nonzeros(pInEachBin(:,10:60,:,i,:));

end

%% Plotting Average Y velocity

plot(1:nt,VelavgY,'-ro')
xlabel('Timestep')
ylabel('average Y-velocity (m/s)')

%(nx,ny,nz,nt,maxP)



%% Plotting a cube
figure
hold on
axis equal

DomainCenter = [.75,.5,.065] ;   % you center point 
Domaindim = [1.5,1,.13] ;  % your cube dimensions 
O = DomainCenter-Domaindim/2 ;       % Get the origin of cube so that P is at center 
plotcube(Domaindim,O,.1,[1 1 1]);   % use function plotcube 
plot3(DomainCenter(1),DomainCenter(2),DomainCenter(3),'*k')
WallCenter = [.5,.55,.065];
Walldim = [.001,.9,.13];
O2 = WallCenter-Walldim/2;
plotcube(Walldim,O2,1,[0 0 0]);
plot3(WallCenter(1),WallCenter(2),WallCenter(3),'*k')

%-----------------------------
u = 1.5/nx
v = 1/ny
w = .13/nz
Cenx = (((InspectXhi-(InspectXlo-1))*u)/2) + (InspectXlo-1)*u ;
Ceny = (((InspectYhi-(InspectYlo-1))*v)/2) + (InspectYlo-1)*v ;
Cenz = (((InspectZhi-(InspectZlo-1))*w)/2) + (InspectZlo-1)*w ;
InspectCenter = [Cenx,Ceny,Cenz];
Inspectdim = [((InspectXhi-(InspectXlo-1))*u), ((InspectYhi-(InspectYlo-1))*v) ((InspectZhi-(InspectZlo-1))*w)];
O3 = InspectCenter-Inspectdim/2
plotcube(Inspectdim,O3,0.6,[0 1 0]);
plot3(InspectCenter(1),InspectCenter(2),InspectCenter(3),'*k')
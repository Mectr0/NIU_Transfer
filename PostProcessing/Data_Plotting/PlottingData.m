%% Plotting Data Script
close all 
clear all
clc
load('PelletData.mat')
fprintf('Number of bins in X: %d  Bin size (m): %f \n', nx, dx)
fprintf('Number of bins in Y: %d  Bin size (m): %f \n', ny, dy)
fprintf('Number of bins in Z: %d  Bin size (m): %f \n', nz, dz)

InspectXlo = 1;     %INPUTS: These Values are inputs for what region the user wants to inspect
InspectXhi = 19;
InspectYlo = 10;
InspectYhi = 61;
InspectZlo = 1;
InspectZhi = 10;
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
u = 1.5/nx;
v = 1/ny;
w = .13/nz;
Cenx = (((InspectXhi-(InspectXlo-1))*u)/2) + (InspectXlo-1)*u ;
Ceny = (((InspectYhi-(InspectYlo-1))*v)/2) + (InspectYlo-1)*v ;
Cenz = (((InspectZhi-(InspectZlo-1))*w)/2) + (InspectZlo-1)*w ;
InspectCenter = [Cenx,Ceny,Cenz];
Inspectdim = [((InspectXhi-(InspectXlo-1))*u), ((InspectYhi-(InspectYlo-1))*v) ((InspectZhi-(InspectZlo-1))*w)];
O3 = InspectCenter-Inspectdim/2;
plotcube(Inspectdim,O3,0.6,[0 1 0]);
plot3(InspectCenter(1),InspectCenter(2),InspectCenter(3),'*k')


%% Plotting Data
for i = 1:nt

Inspect = nonzeros(pInEachBin(:,10:60,:,i,:));
VelavgY(i) = mean(Vel(Inspect,3,i));
Inspect2 = nonzeros(pInEachBin(:,10:60,:,i,:));

end

%% Plotting Average Y velocity
figure
plot(1:nt,VelavgY,'-ro')
xlabel('Timestep')
ylabel('average Y-velocity (m/s)')

%(nx,ny,nz,nt,maxP)


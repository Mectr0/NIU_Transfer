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
InspectYlo = 1;
InspectYhi = 60;
InspectZlo = 1;
InspectZhi = 10;
%% Plotting a cube
%Unless dimensions of simulation change, this shouldn't be touched
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


%% Plotting Quaternion Data
figure
for i = 1:nt
Inspect = nonzeros(pInEachBin(InspectXlo:InspectXhi,InspectYlo:InspectYhi,InspectZlo:InspectZhi,i,:));
%nonzeros is used because the 5th dimension of pInEachBin is preallocated
%with zeros. We don't want these values. 
Q0 = Quat(Inspect,2,i);
Q1 = Quat(Inspect,3,i);
Q2 = Quat(Inspect,4,i);
Q3 = Quat(Inspect,5,i);
histogram(Q1); %Change which quaternion you want
title('Timestep',i)
xlabel('Quaternion q0')
ylabel('Number of Pellets')
pause(0.5)
end

%% Plotting Average Y velocity
for i = 1:nt
Inspect = nonzeros(pInEachBin(InspectXlo:InspectXhi,InspectYlo:InspectYhi,InspectZlo:InspectZhi,i,:));
%nonzeros is used because the 5th dimension of pInEachBin is preallocated
%with zeros. We don't want these values. 
VelavgY(i) = mean(Vel(Inspect,3,i));
end
figure
plot(1:nt,VelavgY,'-ro')
xlabel('Timestep')
ylabel('average Y-velocity (m/s)')

%(nx,ny,nz,nt,maxP)


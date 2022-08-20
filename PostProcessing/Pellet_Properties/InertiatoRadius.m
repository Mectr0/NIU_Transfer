%% Moment of Inertia Script
%This script was used to figure out the size of a simulated pellet made out
%of spheres to have the same moment of inertia of a sample, physical pellet
clear all
close all
clc
%  _  ____     z
%  | |    |    | 
%  | |    |    |
%  h |    |  y o-----x
%  | |    |
%  | |    |
%  -  ----
%
% Iz = (1/2)mr^2
% Ix = Iy = (1/12)m(3r^2 + h^2)
% Isphere = (2/5)mr^2

%% Actual Pellet Dimensions

d = 0.006562 %m  measured from a physical pellet
h = 3.81*d; %roughly 25 mm
rho = 1398.82 % kg/m^3 estimated using caliper and scale

Iz = 0.5*(rho)*(h*pi*(d/2)^2)*(d/2)^2; %Moment of inertia for the real pellet
Ix = (1/12)*(rho)*(h*pi*(d/2)^2)*(3*(d/2)^2 + h^2); %Ix = Iy

%%
R = d/2;
r = sym('r');
m = (rho)*((4/3)*pi*r^3)%Should I keep density the same between sphere and cylinder?
Io = (2/5)*m*r^2; %Moment of inertia for a sphere
Angle = (2*pi)/6; %Six spheres in the outer circle, then one in the middle


%% Iz equivalent
Isum = 0;
for j = 1:10 %For the Iz, each outer sphere is the same distance from the z-axis
    for i = 1:7
        if i == 7 %The seventh one is in the middle
            Isum = Isum + Io;
        else
            Isum = Isum + Io + m*(r)^2;
        end
    end
end
rz= vpasolve(Isum == Iz, r) %,'Real',true)
%0.00165 %m

%% Ix = Iy equivalent
Isum = 0;
for j = 1:5 %I'm using 70 spheres per pellet, using symmetry about x axis
    for i = 1:7
        if i == 7 %Center sphere for each layer
            if j == 1 %The first layer of spheres is only one radius away from the x-axis
            Isum = Isum + Io + m*(r)^2;
            else
            Isum = Isum + Io + m*((1+2*j)*r)^2;
            end
        else
            if j == 1 
            Isum = Isum + Io + m*(sqrt((r)^2 + (r*sin(i*Angle))^2))^2;
            else
            Isum = Isum + Io + m*(sqrt(((1+2*j)*r)^2 + (r*sin(i*Angle))^2))^2; 
            end
        end
    end
end

Isum = 2*Isum; %Symmetry about x-axis
        
rx = vpasolve(Isum == Ix, r) %'Real',true)        

       
        
        
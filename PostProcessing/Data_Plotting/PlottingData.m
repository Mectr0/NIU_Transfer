%% Plotting Data Script
close all 
clear all
clc
load('PelletData1_5.mat')
fprintf('Number of bins in X: %d  Bin size (m): %f \n', nx, dx) %NOTE: These prints the amount of bins in each dimension
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
x = [1; 0; 0];
y = [0; 1; 0];
z = [0; 0; 1]; %Orientation of pellet template
ConveyVec = [1 ; 0 ; 0];
PlotFigTimeStep = 60;
for i = 1:nt
Inspect = nonzeros(pInEachBin(InspectXlo:InspectXhi,InspectYlo:InspectYhi,InspectZlo:InspectZhi,i,:));
%nonzeros is used because the 5th dimension of pInEachBin is preallocated
%with zeros. We don't want these values. 
Q0 = Quat(Inspect,2,i);
Q1 = Quat(Inspect,3,i);
Q2 = Quat(Inspect,4,i);
Q3 = Quat(Inspect,5,i);
OrienVecx = zeros(length(Inspect),3);
OrienVecy = zeros(length(Inspect),3);
OrienVecz = zeros(length(Inspect),3);
Alpha = zeros(length(Inspect),1);
Beta = zeros(length(Inspect),1);
Gamma = zeros(length(Inspect),1);
for j = 1:length(Inspect)
    q0 = Q0(j);
    q1 = Q1(j);
    q2 = Q2(j);
    q3 = Q3(j);
    Ansx = [x(1)*(q0^2+q1^2-q2^2-q3^2)+2*x(2)*(q1*q2-q0*q3)+2*x(3)*((q0*q2)+(q1*q3));   %Finding new orientation for each pellets local axes
           2*x(1)*(q0*q3+q1*q2) + x(2)*(q0^2-q1^2+q2^2-q3^2)+ 2*x(3)*((q2*q3)-(q0*q1));
           2*x(1)*(q1*q3-q0*q2) + 2*x(2)*(q0*q1 + q2*q3) + x(3)*(q0^2-q1^2-q2^2+q3^2)];
    Ansy = [y(1)*(q0^2+q1^2-q2^2-q3^2)+2*y(2)*(q1*q2-q0*q3)+2*y(3)*((q0*q2)+(q1*q3));
           2*y(1)*(q0*q3+q1*q2) + y(2)*(q0^2-q1^2+q2^2-q3^2)+ 2*y(3)*((q2*q3)-(q0*q1));
           2*y(1)*(q1*q3-q0*q2) + 2*y(2)*(q0*q1 + q2*q3) + y(3)*(q0^2-q1^2-q2^2+q3^2)];
    Ansz = [z(1)*(q0^2+q1^2-q2^2-q3^2)+2*z(2)*(q1*q2-q0*q3)+2*z(3)*((q0*q2)+(q1*q3));
           2*z(1)*(q0*q3+q1*q2) + z(2)*(q0^2-q1^2+q2^2-q3^2)+ 2*z(3)*((q2*q3)-(q0*q1));
           2*z(1)*(q1*q3-q0*q2) + 2*z(2)*(q0*q1 + q2*q3) + z(3)*(q0^2-q1^2-q2^2+q3^2)];
    OrienVecx(j,1)= Ansx(1,1);
    OrienVecx(j,2)= Ansx(2,1);
    OrienVecx(j,3)= Ansx(3,1);
    
    OrienVecy(j,1)= Ansy(1,1);
    OrienVecy(j,2)= Ansy(2,1);
    OrienVecy(j,3)= Ansy(3,1);
    
    OrienVecz(j,1)= Ansz(1,1);
    OrienVecz(j,2)= Ansz(2,1);
    OrienVecz(j,3)= Ansz(3,1);
    
    Alpha(j,1) = acosd((ConveyVec(1)*OrienVecx(j,1)+ ConveyVec(2)*OrienVecx(j,2) + ConveyVec(3)*OrienVecx(j,3)));
        
    Beta(j,1) = acosd((ConveyVec(1)*OrienVecy(j,1)+ ConveyVec(2)*OrienVecy(j,2) + ConveyVec(3)*OrienVecy(j,3)));
        
    Gamma(j,1) = acosd((ConveyVec(1)*OrienVecz(j,1)+ ConveyVec(2)*OrienVecz(j,2) + ConveyVec(3)*OrienVecz(j,3)));
    
end
histogram(real(Gamma));
% histogram(Q1); %Change which quaternion you want
% title('Timestep',i)
% xlabel('Quaternion q0')
% ylabel('Number of Pellets')
pause(0.5)
if i == PlotFigTimeStep
saveas(gca,'TestingEPS','epsc')
end
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

%% Plotting Flow Rate
close all
load('../Sorting_Pellets/PelletRate1_5.mat')
load('../Sorting_Pellets/PelletRate2_0.mat')
load('../Sorting_Pellets/PelletRate3_0.mat')
load('../Sorting_Pellets/PelletRate5_0.mat')
load('../Sorting_Pellets/PelletRate8_0.mat')
load('../Sorting_Pellets/PelletRate13_0.mat')
Time = [0:timestep*dumpfreq:timestep*dumpfreq*(length(FlowR)-1)];
% for i = 1:length(Time)
%     if i == 1
%         FlowRate(i) = 0;
%     else
%         FlowRate(i) = (FlowR(i*4))-FlowR(4*(i-1));%/(timestep*dumpfreq);
%     end
% end
figure 
hold on
plot(Time,FlowR1_5,'b*')
plot(Time,FlowR2_0,'gs')
plot(Time,FlowR3_0,'ko')
plot(Time,FlowR5_0,'r*')
plot(Time,FlowR8_0,'bs')
plot(Time,FlowR13_0,'ms')
xlabel('Time (s)')
ylabel('Pellet Number Across Boundary') %Flow Rate (Pellets/Second)')
%legend('1.5 D','2 D','3 D','5 D');
P1_5 = polyfit(Time(60:760), FlowR1_5(60:760), 1);
P2_0 = polyfit(Time(60:760), FlowR2_0(60:760), 1);
P3_0 = polyfit(Time(60:760), FlowR3_0(60:760), 1);
P5_0 = polyfit(Time(60:760), FlowR5_0(60:760), 1);
P8_0 = polyfit(Time(60:760), FlowR8_0(60:760), 1);
P13_0 = polyfit(Time(60:760), FlowR13_0(60:760), 1);
% P1_5 = polyfit(Time, FlowR1_5, 1);
% P2_0 = polyfit(Time, FlowR2_0, 1);
% P3_0 = polyfit(Time, FlowR3_0, 1);
% P5_0 = polyfit(Time, FlowR5_0, 1);

x1_5 = @(x) P1_5(1)*x + P1_5(2);
x2_0 = @(x) P2_0(1)*x + P2_0(2);
x3_0 = @(x) P3_0(1)*x + P3_0(2);
x5_0 = @(x) P5_0(1)*x + P5_0(2);
x8_0 = @(x) P8_0(1)*x + P8_0(2);
x13_0 = @(x) P13_0(1)*x + P13_0(2);
% figure
% hold on
plot(Time, x1_5(Time),'LineWidth',2)
plot(Time, x2_0(Time),'LineWidth',2)
plot(Time, x3_0(Time),'LineWidth',2)
plot(Time, x5_0(Time),'LineWidth',2)
plot(Time, x8_0(Time),'LineWidth',2)
plot(Time, x13_0(Time),'LineWidth',2)


% axis([0 40 -50 200])
% pbaspect([4 2 1]);
fprintf('1.5 D: %f*x + %f\n', P1_5(1), P1_5(2));
fprintf('2.0 D: %f*x + %f\n', P2_0(1), P2_0(2));
fprintf('3.0 D: %f*x + %f\n', P3_0(1), P3_0(2));
fprintf('5.0 D: %f*x + %f\n', P5_0(1), P5_0(2));
fprintf('8.0 D: %f*x + %f\n', P8_0(1), P8_0(2));
fprintf('13.0 D: %f*x + %f\n', P13_0(1), P13_0(2));

%% Plotting Variation

for i = 1:length(FlowR13_0)
    FlowVar13(i) = FlowR13_0(i)-x13_0(Time(i));
    FlowVar8(i) = FlowR8_0(i)-x8_0(Time(i));
    FlowVar5(i) = FlowR5_0(i)-x5_0(Time(i));
    FlowVar3(i) = FlowR3_0(i)-x3_0(Time(i));
    FlowVar2(i) = FlowR2_0(i)-x2_0(Time(i));
    FlowVar1_5(i) = FlowR1_5(i)-x1_5(Time(i));
end
figure 
hold on
% plot(Time(60:end),FlowVar1_5(60:760),'-r*')
%plot(Time(60:end),FlowVar2(60:760),'-bo')
plot(Time(60:end),FlowVar3(60:760),'-ks')
% plot(Time(60:end),FlowVar5(60:760),'-c*')
% plot(Time(60:end),FlowVar8(60:760),'-y*')
% plot(Time(60:end),FlowVar13(60:760),'-g*')
xlabel('Time')
ylabel('Pellet Value - Curve Fit Value')
title('2 Diameter Exit Height')
%plot(Time,FlowVar13,'-r*')


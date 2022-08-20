%Plotting Script created by Connor Murphy through Northern Illinois
%University. 
%     This script was used for collecting experimental results for granular
%     flow experiments under Dr. Nicholas Pohlman and Dr. Jifu Tan.
%     Please follow documentation in each subsection to use this script


%% 1) Setup Step
%*********************************
% Section 1 loads data necessary for ALL following sections
% The main job for this section is loading in data stored in a .mat file. This 
% .mat file needs to be located in the same directory as this script file or
% must include the correct path in the "load" function. The .mat file for this
% section can be created after using 

% /home/mectro/NIU_Undergrad_Research/PostProcessing/Sorting_Pellets

% After running this script, run the second to the last command (currently commented out)
%This should create the desired .mat:

% save('../Data_Plotting/PelletData8_0F.mat','-v7.3')

% NOTE: the "save" command saves all data in the Matlab workspace to a .mat file.
% Make sure your workspace is clear before running the Sorting_Pellets script.
%*********************************

close all 
clear all
clc
load('PelletData5_0F.mat') %Loading the .mat file: contains sorted data and other important properties

fprintf('Number of bins in X: %d  Bin size (m): %f \n', nx, dx) %NOTE: These print the amount of bins in each dimension
fprintf('Number of bins in Y: %d  Bin size (m): %f \n', ny, dy)
fprintf('Number of bins in Z: %d  Bin size (m): %f \n', nz, dz)

PelletMass =0.000294; %kg
InspectXlo = 1;     %INPUTS: These Values are inputs for what region the user wants to inspect
InspectXhi = nx;    %nx, ny , and nz are the number of bins in each direction
InspectYlo = 1;     %The .mat file should have dx, dy, and dz values. These are the dimensions of each bin in (m)
InspectYhi = ny;    %Use this to figure out what part of the 0.5 m (X) by 1 m (Y) by 0.13 m (z) domain you want to inspect 
InspectZlo = 1;
InspectZhi = nz;
%% 2) Plotting a cube
%Unless dimensions of simulation change, this shouldn't be touched
%This is code that should show the user what part of the domain box he will
%be inspecting
%NOTE: Section 2 needs section 1 to be run first

%REQUIREMENTS: Needs Section 1 to be run beforehand

figure
hold on
axis equal

DomainCenter = [.75,.5,.065] ;   % your center point 
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
u = 0.5/nx;
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

%% 3) Quaternion calculations 
%************************************
% Section 3 uses sorted quaternion data aquired from section 1 and uses it to
% generate contour plots that show the average angle between pellets and the 
% motion of the conveyor belt as a function of time and X position
%NOTE: Section 3 needs section 1 to be run first
%************************************

%REQUIREMENTS: Needs Section 1 to be run beforehand

ConveyVec = [1 ; 0 ; 0]; %Conveyor moves in the X-direction
x = [1; 0; 0];
y = [0; 1; 0];
z = [0; 0; 1]; %Orientation of pellet template
GammaMean = zeros(InspectXhi-1,1);
it = 10; %The timestep we want to look at
iy = 7; %Highest bin in y direction we want to look at.

%************************************
%SUBJECT TO USER CHANGE
%NOTE:(dy = 0.0195m)
%Height to top of flight: 0.0352m  Inspecting two bins (1:2): 0.039 m  
%Inspecting above flights (3:10) : Region(0.039 to 0.195 m)
%Inspecting top flights (11:22) : Region(0.195 to 0.429m) 
%************************************
for i = 1:nt
    for j = 1:InspectXhi-1
        Inspect = nonzeros(pInEachBin(j,5:6,1:nz,i,:));%NOTE: Explanation of pInEachBin is found in SortingScript.m
        %nonzeros is used because the 5th dimension of pInEachBin is preallocated
        %with zeros. We don't want these values.
        Gamma = zeros(length(Inspect),1);
        for k = 1:length(Inspect)
            idx = find(Quat(:,1,i)== Inspect(k));
            q0 = Quat(idx,2,i);
            q1 = Quat(idx,3,i);
            q2 = Quat(idx,4,i);
            q3 = Quat(idx,5,i);
            Ansz = [z(1)*(q0^2+q1^2-q2^2-q3^2)+2*z(2)*(q1*q2-q0*q3)+2*z(3)*((q0*q2)+(q1*q3));
                    2*z(1)*(q0*q3+q1*q2) + z(2)*(q0^2-q1^2+q2^2-q3^2)+ 2*z(3)*((q2*q3)-(q0*q1));
                    2*z(1)*(q1*q3-q0*q2) + 2*z(2)*(q0*q1 + q2*q3) + z(3)*(q0^2-q1^2-q2^2+q3^2)];
            OrienVecz(k,1)= Ansz(1,1);
            OrienVecz(k,2)= Ansz(2,1);
            OrienVecz(k,3)= Ansz(3,1);
            Gamma(k,1) = acosd((ConveyVec(1)*OrienVecz(k,1)+ ConveyVec(2)*OrienVecz(k,2) + ConveyVec(3)*OrienVecz(k,3)));
            if Gamma(k,1) < 0
                fprint('NEGATIVE ANGLE')
            end
            if Gamma(k,1) > 90
                Gamma(k,1) = 180 - Gamma(k,1);
            end
            if Inspect(k) ~= idx
                fprintf('WARNING: Pellet ID does not match row number (Quat Calculation p2)')
                %Just in case the sorted rows in Quat dont have id values that
                %match up to the row number. Shouldn't be a problem, just
                %checking.
            end

        end

        Gammastd(j,i) = std(Gamma(:,1));
        GammaMean(j,i) = mean(Gamma(:,1));
    end
end 
close all
figure
colormap('jet')
%surf(real(GammaMean));
Time = [0:timestep*dumpfreq:timestep*dumpfreq*(length(FlowR)-1)]; %Multiply timestep with dump frequency as well as the total
%number of timesteps will give you how many seconds the simulation
%represents
Pos = [0:0.0263:0.5]; %This is subject to change if the user changes the domain (0.5m/19) = 0.0263 Needed to reduce from 21 points to 20
contourf(Time,Pos,real(GammaMean(1:InspectXhi-1,:)), [0:2:90],'edgecolor','none')
caxis([22 80]) %This readjusts the colorbar. Make sure no values are outside this range
colorbar
shading INTERP

%------------------Use this to overlay belt speed on contour
hold on
BeltSpeed = 0.02653; %m/s  

for ii=1:19
    BeltTime = [0*timestep*dumpfreq, (66*(ii-1)+78)*timestep*dumpfreq];
    BeltPos = @(n) BeltSpeed*n + 0.478535-(ii-1)*0.0813;
    plot(BeltTime,BeltPos(BeltTime), '--r', 'LineWidth',2);
end
% annotation('doublearrow',[0.288 0.4], [0.925 0.925]) %0.132015 0.219645
% annotation('doublearrow',[0.4 0.52], [0.925 0.925])
% annotation('doublearrow',[0.52 0.63], [0.925 0.925])
% BeltPitch=0.08763;
% FlightPos1=[0:0.0001*500*20*BeltSpeed:0.5];
% FlightPos2=[0+BeltPitch:timestep*dumpfreq*20*BeltSpeed:0.5];
% FlightTime=[5:23];
% plot(FlightTime,FlightPos1,'bs','MarkerFaceColor','b')
% plot(FlightTime(1:16),FlightPos2,'ro','MarkerFaceColor','r')
%------------------

xlabel('Time (s)','Fontsize',16), ylabel('X Position (m)','Fontsize',16);
view(90, 90);
 axis([0 37.95 0 0.5])
%saveas(gca,'Contour5D5to6.eps','epsc')
min(min(real(GammaMean(1:nx,:))))%checking to see if my colorbar range is good
max(max(real(GammaMean(1:nx,:))))


%% 4) Plotting Average Y velocity

%REQUIREMENTS: Needs Section 1 to be run beforehand

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

%% 5) Plotting Flow Rate

%REQUIREMENTS: Needs Section 1 to be run beforehand

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
box on

yyaxis left
ylabel('Expelled Pellets') %Flow Rate (Pellets/Second)')
xlabel('Time (s)')
axis([0 40 -100 10000])
%subplot(2,2,4),plot(Time,FlowR13_0,'ms')
plot(Time,FlowR13_0,'--','LineWidth', 2, 'color','#FF00FF')
P13_0 = polyfit(Time(60:760), FlowR13_0(60:760), 1);
x13_0 = @(x) P13_0(1)*x + P13_0(2);
% plot(Time, x13_0(Time),'LineWidth',2)

%subplot(2,2,4),plot(Time,FlowR8_0,'bs')
plot(Time,FlowR8_0,'-.','LineWidth', 2, 'color','#0000FF')
P8_0 = polyfit(Time(60:760), FlowR8_0(60:760), 1);
x8_0 = @(x) P8_0(1)*x + P8_0(2);
%plot(Time, x8_0(Time),'LineWidth',2)

%subplot(2,2,4),plot(Time,FlowR5_0,'r*')
plot(Time,FlowR5_0,':','LineWidth', 2, 'color','#77AC30')
P5_0 = polyfit(Time(60:760), FlowR5_0(60:760), 1);
x5_0 = @(x) P5_0(1)*x + P5_0(2);
% plot(Time, x5_0(Time),'LineWidth',2)

%subplot(2,2,3),plot(Time,FlowR3_0,'ko')
plot(Time,FlowR3_0,'-.','LineWidth', 2, 'color','#0072BD')
P3_0 = polyfit(Time(60:760), FlowR3_0(60:760), 1);
x3_0 = @(x) P3_0(1)*x + P3_0(2);
% plot(Time, x3_0(Time),'LineWidth',2)

%subplot(2,2,2),plot(Time,FlowR2_0,'gs')
plot(Time,FlowR2_0,'-','LineWidth', 2, 'color','#000000')
P2_0 = polyfit(Time(60:760), FlowR2_0(60:760), 1);
x2_0 = @(x) P2_0(1)*x + P2_0(2);
% plot(Time, x2_0(Time),'LineWidth',2)

%subplot(2,2,1), plot(Time,FlowR1_5,'b*')
plot(Time,FlowR1_5,'-','LineWidth', 2, 'color','#A2142F')
P1_5 = polyfit(Time(60:760), FlowR1_5(60:760), 1);
x1_5 = @(x) P1_5(1)*x + P1_5(2);
% plot(Time, x1_5(Time),'LineWidth',2)

ylabel('Expelled Pellets','Fontsize',16) %Flow Rate (Pellets/Second)')
xlabel('Time (s)','Fontsize',16)
leg = legend('13 D','8 D','5 D','3 D','2 D','1.5 D');
title(leg,'Exit Height in Pellet Diameter (D)')
%axis([0 40 -100 8600])
%axis([0 7 -100 1450])

yyaxis right
%axis([0 40 -100*PelletMass 8600*PelletMass])
%axis([0 7 -100*PelletMass 1450*PelletMass])
ylim([0 , 16.5])
Mass_flow_Data;
dt = 2;
time8D = [0:dt:dt*(length(data8D)-1)]';
time5D = [0:dt:dt*(length(data5D)-1)]';
time3D = [0:dt:dt*(length(data3D)-1)]';
time2D = [0:dt:dt*(length(data2D)-1)]';
plot(time8D,data8D,'^','Markersize',6, 'color','#0000FF','MarkerFaceColor','#0000FF');
plot(time5D,data5D,'>','Markersize',6, 'color','#77AC30','MarkerFaceColor','#77AC30');
plot(time3D,data3D,'rv','Markersize',6,'Color','#0072BD','MarkerFaceColor','#0072BD');
plot(time2D,data2D,'<','Markersize',6,'color','#000000','MarkerFaceColor','#000000');
ylabel('Mass Expelled (kg)','Fontsize',16)
leg = legend('13 D','8 D','5 D','3 D','2 D','1.5 D','Exp. $m$','Location','north','Orientation','Horizontal','NumColumns',4,'Interpreter','latex');
legend boxoff
title(leg,'Exit Height in Pellet Diameter (D)')


% p = get(gca, 'Position');
% h = axes('Parent', gcf, 'Position', [p(1)+.06 p(2)+.06 p(3)-.5 p(4)-.5]);
% plot(h,Time,FlowR13_0,'--','LineWidth', 2, 'color','#FF00FF')

axes('Position',[.2 .5 .25 .25])
box on
hold on
plot(Time(1:60),FlowR13_0(1:60),'--','LineWidth', 2, 'color','#FF00FF')
plot(Time(1:60),FlowR8_0(1:60),'-.','LineWidth', 2, 'color','#0000FF')
plot(Time(1:60),FlowR5_0(1:60),':','LineWidth', 2, 'color','#77AC30')
plot(Time(1:60),FlowR3_0(1:60),'-.','LineWidth', 2, 'color','#0072BD')
plot(Time(1:60),FlowR2_0(1:60),'-','LineWidth', 2, 'color','#000000')
plot(Time(1:60),FlowR1_5(1:60),'-','LineWidth', 2, 'color','#A2142F')

FlowRates = [P1_5(1) P2_0(1) P3_0(1) P5_0(1) P8_0(1) P13_0(1)];

%saveas(gca,'FlowRateCombine.eps','epsc')

% pbaspect([4 2 1]);
fprintf('1.5 D: %f*x + %f\n', P1_5(1), P1_5(2));
fprintf('2.0 D: %f*x + %f\n', P2_0(1), P2_0(2));
fprintf('3.0 D: %f*x + %f\n', P3_0(1), P3_0(2));
fprintf('5.0 D: %f*x + %f\n', P5_0(1), P5_0(2));
fprintf('8.0 D: %f*x + %f\n', P8_0(1), P8_0(2));
fprintf('13.0 D: %f*x + %f\n', P13_0(1), P13_0(2));

%% 6) Plotting Variation
close all

%REQUIREMENTS: Needs Sections 1 and 5 to be run beforehand

for i = 1:length(FlowR13_0)
    FlowVar13(i) = FlowR13_0(i)-x13_0(Time(i));
    FlowVar8(i) = FlowR8_0(i)-x8_0(Time(i));
    FlowVar5(i) = FlowR5_0(i)-x5_0(Time(i));
    FlowVar3(i) = FlowR3_0(i)-x3_0(Time(i));
    FlowVar2(i) = FlowR2_0(i)-x2_0(Time(i));
    FlowVar1_5(i) = FlowR1_5(i)-x1_5(Time(i));
end
figure 
box on
hold on
plot(Time(60:end),FlowVar13(60:760),'--','LineWidth', 2, 'color','#FF00FF')
plot(Time(60:end),FlowVar8(60:760),'-.','LineWidth', 2, 'color','#0000FF')
plot(Time(60:end),FlowVar5(60:760),':','LineWidth', 2, 'color','#77AC30')
plot(Time(60:end),FlowVar3(60:760),'-.','LineWidth', 2, 'color','#0072BD')
plot(Time(60:end),FlowVar2(60:760),'-','LineWidth', 2, 'color','#000000')
plot(Time(60:end),FlowVar1_5(60:760),'-','LineWidth', 2, 'color','#A2142F')
yline(0,'-r', 'LineWidth', 3)
xlabel('Time (s)','Fontsize',16)
ylabel('Pellet Value - Curve Fit Value','Fontsize',16)
%title('2 Diameter Exit Height')
leg = legend('13 D','8 D','5 D','3 D','2 D','1.5 D','Location','north','Orientation','Horizontal','NumColumns',4)
title(leg,'Exit Height in Pellet Diameter (D)')
legend Boxoff
axis([0 40 -110 60])
saveas(gca,'FlowRVariation.eps','epsc')

figure

t_end = 37.95-2.95; %2.95 because the first value is 0 (instead of using a value of 3)
df = 1/t_end; %frequency increment
fs = 1/(timestep*dumpfreq); %Sampling frequency
ft5 = fft(FlowVar5(61:760));
ft5abs = sqrt(imag(ft5).^2 + real(ft5).^2)/(760-60);
f = linspace(0,fs,fs/df);
STEM = stem(f,ft5abs,'LineWidth', 2,'color','#FF00FF');
set(gca,'fontsize',16)
xlabel('Frequency (Hz)','Fontsize',16)
ylabel('Pellet Flow Variance','Fontsize',16)

hold on
for i = 1:9
xline(0.3264*i,'--r','linewidth', 1.5);
end
xlim([-0.1 , 3])
ylim([0 , 14])
%breakyaxis([6 13])
%saveas(gca,'STD.eps','epsc')
%% 7) Standard Deviation, Skewness, and Kurtosis
close all

%REQUIREMENTS: Needs Sections 1, 5 and 6 to be run beforehand

% Sigma1_5 = std(FlowVar1_5(60:760)); %Sample std (normalized by n-1)
% Sigma2 = std(FlowVar2(60:760));
% Sigma3 = std(FlowVar3(60:760));
% Sigma5 = std(FlowVar5(60:760));
% Sigma8 = std(FlowVar8(60:760));
% Sigma13 = std(FlowVar13(60:760));
Sigma(1)=std(FlowVar1_5(60:760)); 
Sigma(2)=std(FlowVar2(60:760));
Sigma(3)=std(FlowVar3(60:760));
Sigma(4)=std(FlowVar5(60:760));
Sigma(5)=std(FlowVar8(60:760));
Sigma(6)=std(FlowVar13(60:760)); 

% Skew1_5 = skewness(FlowVar1_5(60:760)); %Should we correct bias? Ex: skewness(X,0)
% Skew2 = skewness(FlowVar2(60:760));
% Skew3 = skewness(FlowVar3(60:760));
% Skew5 = skewness(FlowVar5(60:760));
% Skew8 = skewness(FlowVar8(60:760));
% Skew13 = skewness(FlowVar13(60:760
Skew(1)=skewness(FlowVar1_5(60:760));
Skew(2)=skewness(FlowVar2(60:760));
Skew(3)=skewness(FlowVar3(60:760));
Skew(4)=skewness(FlowVar5(60:760));
Skew(5)=skewness(FlowVar8(60:760));
Skew(6)=skewness(FlowVar13(60:760));

% Kurt1_5 = kurtosis(FlowVar1_5(60:760));%Should we correct bias? Ex: kurtosis(X,0)
% Kurt2 = kurtosis(FlowVar2(60:760));
% Kurt3 = kurtosis(FlowVar3(60:760));
% Kurt5 = kurtosis(FlowVar5(60:760));
% Kurt8 = kurtosis(FlowVar8(60:760));
% Kurt13 = kurtosis(FlowVar13(60:760));
Kurt(1)=kurtosis(FlowVar1_5(60:760));
Kurt(2)=kurtosis(FlowVar2(60:760));
Kurt(3)=kurtosis(FlowVar3(60:760));
Kurt(4)=kurtosis(FlowVar5(60:760));
Kurt(5)=kurtosis(FlowVar8(60:760));
Kurt(6)=kurtosis(FlowVar13(60:760));

SimList = [1.5,2,3,5,8,13];

Error = figure;
%plot(SimList,FlowRates, 'bv','MarkerFaceColor', 'b')
Error = errorbar(SimList-0.1,FlowRates,Sigma,'ko','MarkerFaceColor', 'k','LineWidth', 1.5)

colororder({'k','r'})
grid on
yyaxis left
xlabel('Exit Height (D)','Fontsize',16)
ylabel('Pellet Flow Rate (pellets/s)','Fontsize',16)
set(gca, 'XTick',SimList) 
labelpoints(SimList(3:5),FlowRates(3:5),round(Sigma(3:5)),'E', 0.2, 1, 'Color', 'k','Fontsize', 14)
labelpoints(SimList(1),FlowRates(1),round(Sigma(1)),'N', 0.65, 1, 'Color', 'k','Fontsize', 14)
labelpoints(SimList(2),FlowRates(2),round(Sigma(2)),'N', 0.9, 1, 'Color', 'k','Fontsize', 14)
labelpoints(SimList(6),FlowRates(6),round(Sigma(6)),'W', 0.2, 1, 'Color', 'k','Fontsize', 14)
axis([1 13.5 0 275]);

yyaxis right

Xexp = [2.1, 3.1, 5.1, 8.1];
Yexp = [0.110968, 0.125113, 0.156776, 0.273831];
errorexp = [0.046, 0.038, 0.058, 0.20]
errorbar(Xexp,Yexp,errorexp,'rv','MarkerFaceColor', 'r','LineWidth', 1.5)
labelpoints(Xexp(2),Yexp(2),errorexp(2),'N', 1.3, 1, 'Color', 'r','Fontsize', 14)
labelpoints(Xexp(3),Yexp(3),errorexp(3),'N', 1.7, 1, 'Color', 'r','Fontsize', 14)
labelpoints(Xexp(1),Yexp(1),errorexp(1),'S', 1, 1, 'Color', 'r','Fontsize', 14)
labelpoints(Xexp(4),Yexp(4),errorexp(4),'S', 4, 1, 'Color', 'r','Fontsize', 14)

axis([1 13.5 0 0.5]);
ylabel('Mass Flow Rate (kg/s)','Fontsize',16)
% figure 
% plot(SimList,Skew, ':ro','MarkerFaceColor', 'r')
% grid on
% xlabel('Exit Height (D)','Fontsize',14)
% ylabel('Skewness ','Fontsize',14)
% 
% figure
% plot(SimList,Kurt, ':ks','MarkerFaceColor', 'k')
% grid on
% xlabel('Exit Height (D)','Fontsize',14)
% ylabel('Kurtosis ','Fontsize',14)

%% 8) Velocity Quiver Plot
close all

%REQUIREMENTS: Needs Section 1 to be run beforehand

time = 560;
vxmean = zeros(InspectYhi, InspectXhi);
vymean = zeros(InspectYhi, InspectXhi);
for j = 1:InspectYhi
    for i = 1:InspectXhi
        Inspect = nonzeros(pInEachBin(i,j,1:nz,time,:));
        %nonzeros is used because the 5th dimension of pInEachBin is preallocated
        %with zeros. We don't want these values.
        vx = zeros(length(Inspect),1);
        vy = zeros(length(Inspect),1);
        for k = 1:length(Inspect)
            idx = find(Vel(:,1,time)== Inspect(k));
            vx(k,1) = Vel(idx,2,time);
            vy(k,1) = Vel(idx,3,time);
            if Inspect(k) ~= idx
                fprintf('WARNING: Pellet ID does not match row number (Quat Calculation p2)')
                %Just in case the sorted rows in Quat dont have id values that
                %match up to the row number. Shouldn't be a problem, just
                %checking.
            end

        end
        vxmean(j,i) = mean(vx);
        vymean(j,i) = mean(vy);
    end
end
figure
scale = 4;
quiver([dx/2:dx:nx*dx],[dy/2:dy:ny*dy], vxmean*scale, vymean*scale, 'Autoscale', 'off');
hold on
plot([.5,.5], [(0.0098+0.0254+13*0.00738),0.5],'-','LineWidth', 2, 'color','#A2142F');
quiver(0.4,0.03,(.02653)*scale,0*scale,'LineWidth', 4, 'color','#A2142F', 'Autoscale', 'off', 'MaxHeadSize', .75)
xlabel('X Location (m)','Fontsize',16)
ylabel('Y Location (m)','Fontsize',16)
axis equal
axis([0 0.6 0 0.4])
%saveas(gca,'13DQuiver.eps','epsc')


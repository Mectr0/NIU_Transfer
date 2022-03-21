%% Plotting Data Script
close all 
clear all
clc
load('PelletData13_0.mat')
fprintf('Number of bins in X: %d  Bin size (m): %f \n', nx, dx) %NOTE: These prints the amount of bins in each dimension
fprintf('Number of bins in Y: %d  Bin size (m): %f \n', ny, dy)
fprintf('Number of bins in Z: %d  Bin size (m): %f \n', nz, dz)

InspectXlo = 1;     %INPUTS: These Values are inputs for what region the user wants to inspect
InspectXhi = nx;    
InspectYlo = 1;
InspectYhi = ny;
InspectZlo = 1;
InspectZhi = nz;
%% Plotting a cube
%Unless dimensions of simulation change, this shouldn't be touched

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
title('Timestep',i)
xlabel('Angle Gamma (degrees)')
ylabel('Number of Pellets')
pause(0.5)
if i == PlotFigTimeStep
saveas(gca,'TestingEPS','epsc')
end
end
%% Quaternion calculations part 2

ConveyVec = [1 ; 0 ; 0];
x = [1; 0; 0];
y = [0; 1; 0];
z = [0; 0; 1]; %Orientation of pellet template
GammaMean = zeros(InspectXhi,1);
it = 760; %The timestep we are interested in
iy = 10; %Highest bin in y direction we want to look at. (dy*10 = 26 pellet diameters)
for j = 1:InspectXhi
    Inspect = nonzeros(pInEachBin(j,1:iy,1:nz,it,:));
    %nonzeros is used because the 5th dimension of pInEachBin is preallocated
    %with zeros. We don't want these values.
    Gamma = zeros(length(Inspect),1);
    for k = 1:length(Inspect)
        idx = find(Quat(:,1,it)== Inspect(k));
        q0 = Quat(idx,2,it);
        q1 = Quat(idx,3,it);
        q2 = Quat(idx,4,it);
        q3 = Quat(idx,5,it);
        Ansz = [z(1)*(q0^2+q1^2-q2^2-q3^2)+2*z(2)*(q1*q2-q0*q3)+2*z(3)*((q0*q2)+(q1*q3));
                2*z(1)*(q0*q3+q1*q2) + z(2)*(q0^2-q1^2+q2^2-q3^2)+ 2*z(3)*((q2*q3)-(q0*q1));
                2*z(1)*(q1*q3-q0*q2) + 2*z(2)*(q0*q1 + q2*q3) + z(3)*(q0^2-q1^2-q2^2+q3^2)];
        OrienVecz(k,1)= Ansz(1,1);
        OrienVecz(k,2)= Ansz(2,1);
        OrienVecz(k,3)= Ansz(3,1);
        Gamma(k,1) = acosd((ConveyVec(1)*OrienVecz(k,1)+ ConveyVec(2)*OrienVecz(k,2) + ConveyVec(3)*OrienVecz(k,3)));
        if Inspect(k) ~= idx
            fprintf('WARNING: Pellet ID does not match row number (Quat Calculation p2)')
            %Just in case the sorted rows in Quat dont have id values that
            %match up to the row number. Shouldn't be a problem, just
            %checking.
        end

    end
    Gammastd(j) = std(Gamma(:,1));
    GammaMean(j) = mean(Gamma(:,1));
end
figure
Spacing = [0:dx:(dx*nx)-dx];
GammaMeanPlot = GammaMean(:,1);
%plot(Spacing, GammaMeanPlot, 'k-o')
errorbar(Spacing,GammaMeanPlot, Gammastd, 'k-o')
xlabel = ('Discretized position along hopper (m)')
ylabel = ('Angle between length of pellet and belt motion')
title('Timestep:',it)

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
box on

%subplot(2,2,4),plot(Time,FlowR13_0,'ms')
plot(Time,FlowR13_0,'--','LineWidth', 2, 'color','#FF00FF')
P13_0 = polyfit(Time(60:760), FlowR13_0(60:760), 1);
x13_0 = @(x) P13_0(1)*x + P13_0(2);
% plot(Time, x13_0(Time),'LineWidth',2)

%subplot(2,2,4),plot(Time,FlowR8_0,'bs')
plot(Time,FlowR8_0,'-.','LineWidth', 2, 'color','#0000FF')
P8_0 = polyfit(Time(60:760), FlowR8_0(60:760), 1);
x8_0 = @(x) P8_0(1)*x + P8_0(2);
% plot(Time, x8_0(Time),'LineWidth',2)

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
plot(Time,FlowR1_5,'--','LineWidth', 2, 'color','#A2142F')
P1_5 = polyfit(Time(60:760), FlowR1_5(60:760), 1);
x1_5 = @(x) P1_5(1)*x + P1_5(2);
% plot(Time, x1_5(Time),'LineWidth',2)

xlabel('Time (s)')
ylabel('Expelled Pellets') %Flow Rate (Pellets/Second)')
leg = legend('13 D','8 D','5 D','3 D','2 D','1.5 D');
title(leg,'Exit Height in Pellet Diameter (D)')

axis([0 40 -100 6500])

%saveas(gca,'FlowRateCombine','epsc')

% pbaspect([4 2 1]);
fprintf('1.5 D: %f*x + %f\n', P1_5(1), P1_5(2));
fprintf('2.0 D: %f*x + %f\n', P2_0(1), P2_0(2));
fprintf('3.0 D: %f*x + %f\n', P3_0(1), P3_0(2));
fprintf('5.0 D: %f*x + %f\n', P5_0(1), P5_0(2));
fprintf('8.0 D: %f*x + %f\n', P8_0(1), P8_0(2));
fprintf('13.0 D: %f*x + %f\n', P13_0(1), P13_0(2));

%% Plotting Variation
close all
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
plot(Time(60:end),FlowVar13(60:760),'--','LineWidth', 2, 'color','#FF00FF')
plot(Time(60:end),FlowVar8(60:760),'-.','LineWidth', 2, 'color','#0000FF')
plot(Time(60:end),FlowVar5(60:760),':','LineWidth', 2, 'color','#77AC30')
plot(Time(60:end),FlowVar3(60:760),'-.','LineWidth', 2, 'color','#0072BD')
plot(Time(60:end),FlowVar2(60:760),'-','LineWidth', 2, 'color','#000000')
plot(Time(60:end),FlowVar1_5(60:760),'--','LineWidth', 2, 'color','#A2142F')

xlabel('Time (S)')
ylabel('Pellet Value - Curve Fit Value')
%title('2 Diameter Exit Height')
leg = legend('13 D','8 D','5 D','3 D','2 D','1.5 D')
title(leg,'Exit Height in Pellet Diameter (D)')
axis([0 40 -110 60])
%plot(Time,FlowVar13,'-r*')

%Variance Histograms
figure
subplot(3,2,1),histogram(FlowVar13(60:760),40)
title('Variance for 13 D')
subplot(3,2,2),histogram(FlowVar8(60:760),40)
title('Variance for 8 D')
subplot(3,2,3),histogram(FlowVar5(60:760),40)
title('Variance for 5 D')
subplot(3,2,4),histogram(FlowVar3(60:760),40)
title('Variance for 3 D')
subplot(3,2,5),histogram(FlowVar2(60:760),40)
title('Variance for 2 D')
subplot(3,2,6),histogram(FlowVar1_5(60:760),40)
title('Variance for 1.5 D')

%Fourier Transform
figure
t_end = 37.95-2.95; %2.95 because the first value is 0
df = 1/t_end; %frequency increment
fs = 1/(timestep*dumpfreq); %Sampling frequency
ft13 = fft(FlowVar13(61:760));
ft13abs = sqrt(imag(ft13).^2 + real(ft13).^2);
f = linspace(0,fs,fs/df);
stem(f,ft13abs);

%% Standard Deviation, Skewness, and Kurtosis
close all

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

figure
plot(SimList,Sigma, ':bv','MarkerFaceColor', 'b')
grid on
xlabel('Exit Height (D)')
ylabel('Standard Deviation of Flowrate Variance')

figure 
plot(SimList,Skew, ':ro','MarkerFaceColor', 'r')
grid on
xlabel('Exit Height (D)')
ylabel('Skewness of Flowrate Variance Distribution')

figure
plot(SimList,Kurt, ':ks','MarkerFaceColor', 'k')
grid on
xlabel('Exit Height (D)')
ylabel('Kurtosis of Flowrate Variance Distribution')




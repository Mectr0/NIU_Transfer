%% Orientation Quat Post Process
%-----------------------------------------------------

%Connor Murphy

%This code takes a dump file for multisphere that has quaternion
%and center of mass values:

% q0 q1 q2 q3 cmx cmy cmz


%-------------------------------------------------------
clear all
numPellets = 300;
%partsperpellet = 70;
numTimesteps = 50000 ; 
dumpFrequency = 1000 ;
delta = 0.006562 % m   Pellet diameter

DumpNumber = numTimesteps/dumpFrequency;

name = 'dump.pellet_orien';
fid = fopen(name,'r+');
Quat =  zeros(DumpNumber, 4, numPellets);
PosCM = zeros(DumpNumber, 3, numPellets);
for k = 1:5
        fgetl(fid); % The dump files comes with 18 lines of script before data starts
end
 
Xlo = fscanf(fid,'%g',1);
Xhi = fscanf(fid,'%g',1);
Ylo = fscanf(fid,'%g',1);
Yhi = fscanf(fid,'%g',1);
Zlo = fscanf(fid,'%g',1);
Zhi = fscanf(fid,'%g',1);
fgetl(fid);
fgetl(fid);

for i = 1:DumpNumber
    for k = 1:9 %9 lines of script between each dumped timestep
        fgetl(fid);
    end
    for j = 1:numPellets %Reading data from each line 
        Quat(i,1,j)= fscanf(fid,'%g',1);
        Quat(i,2,j)= fscanf(fid,'%g',1);
        Quat(i,3,j)= fscanf(fid,'%g',1);
        Quat(i,4,j)= fscanf(fid,'%g',1);
        PosCM(i,1,j)= fscanf(fid,'%g',1);
        PosCM(i,2,j)= fscanf(fid,'%g',1);
        PosCM(i,3,j)= fscanf(fid,'%g',1);
    end
    fgetl(fid);
    
end
%Orientation vectors for each pellet's local axis. 
%The pellet template is positioned lengthwise along the z axis

nx = floor((Xhi-(-0.5))/(1.7*delta))   %floor((Xhi-Xlo)/(1.7*delta)) The x domain is too long
ny = floor((Yhi-Ylo)/(.7*delta))   
nz = floor((Zhi-Zlo)/(1.7*delta))

%--------------------------------------------
Domainx = [1 134]; %INPUT: REGION OF INTEREST
if Domainx(1,2) > nx || Domainx(1,1) < 1
    error('Region of interest in x must match domain mesh')
end
%--------------------------------------------
Domainy = [1 217]; %INPUT: REGION OF INTEREST
if Domainy(1,2) > ny || Domainy(1,1) < 1
    error('Region of interest in y must match domain mesh')
end
%--------------------------------------------
Domainz = [1 13];  %INPUT: REGION OF INTEREST           
if Domainz(1,2) > nz || Domainz(1,1) < 1
    error('Region of interest in z must match domain mesh')
end

for kk = 1:(Domainz(1,2)-Domainz(1,1))
    for jj = 1:(Domainy(1,2)-Domainy(1,1))
        for ii = 1:(Domainx(1,2)-Domainx(1,1))
        end
    end
end
OrienVecx = zeros(DumpNumber, 3, numPellets); %Allocating space 
OrienVecy = zeros(DumpNumber, 3, numPellets);
OrienVecz = zeros(DumpNumber, 3, numPellets);
x = [1; 0; 0];
y = [0; 1; 0];
z = [0; 0; 1];

for i = 1:DumpNumber
    for j = 1:numPellets
    q0 = Quat(i,1,j); %Using the quaternion for each timestep
    q1 = Quat(i,2,j);
    q2 = Quat(i,3,j);
    q3 = Quat(i,4,j);
 
    Ansx = [x(1)*(q0^2+q1^2-q2^2-q3^2)+2*x(2)*(q1*q2-q0*q3)+2*x(3)*((q0*q2)+(q1*q3));   %Finding new orientation for each pellets local axes
           2*x(1)*(q0*q3+q1*q2) + x(2)*(q0^2-q1^2+q2^2-q3^2)+ 2*x(3)*((q2*q3)-(q0*q1));
           2*x(1)*(q1*q3-q0*q2) + 2*x(2)*(q0*q1 + q2*q3) + x(3)*(q0^2-q1^2-q2^2+q3^2)];
    
    Ansy = [y(1)*(q0^2+q1^2-q2^2-q3^2)+2*y(2)*(q1*q2-q0*q3)+2*y(3)*((q0*q2)+(q1*q3));
           2*y(1)*(q0*q3+q1*q2) + y(2)*(q0^2-q1^2+q2^2-q3^2)+ 2*y(3)*((q2*q3)-(q0*q1));
           2*y(1)*(q1*q3-q0*q2) + 2*y(2)*(q0*q1 + q2*q3) + y(3)*(q0^2-q1^2-q2^2+q3^2)];
    
    Ansz = [z(1)*(q0^2+q1^2-q2^2-q3^2)+2*z(2)*(q1*q2-q0*q3)+2*z(3)*((q0*q2)+(q1*q3));
           2*z(1)*(q0*q3+q1*q2) + z(2)*(q0^2-q1^2+q2^2-q3^2)+ 2*z(3)*((q2*q3)-(q0*q1));
           2*z(1)*(q1*q3-q0*q2) + 2*z(2)*(q0*q1 + q2*q3) + z(3)*(q0^2-q1^2-q2^2+q3^2)];
    
%Storing new orientation vectors    
    OrienVecx(i,1,j)= Ansx(1,1);
    OrienVecx(i,2,j)= Ansx(2,1);
    OrienVecx(i,3,j)= Ansx(3,1);
    
    OrienVecy(i,1,j)= Ansy(1,1);
    OrienVecy(i,2,j)= Ansy(2,1);
    OrienVecy(i,3,j)= Ansy(3,1);
    
    OrienVecz(i,1,j)= Ansz(1,1);
    OrienVecz(i,2,j)= Ansz(2,1);
    OrienVecz(i,3,j)= Ansz(3,1);
    end
end

ConveyVec = [1 ; 0 ; 0];
%Allocating for angles between each orientation vector and the conveyor
%belt (moving in 1i + 0j + 0k)
Alpha = zeros(DumpNumber, 1, numPellets);

Beta = zeros(DumpNumber, 1, numPellets);

Gamma = zeros(DumpNumber, 1, numPellets);
for i = 1:DumpNumber
    for j = 1:numPellets
        %Using the dot product
        %|ConveyVec|*|OrienVec| will always be one
        Alpha(i,1,j) = acosd((ConveyVec(1)*OrienVecx(i,1,j)+ ConveyVec(2)*OrienVecx(i,2,j) + ConveyVec(3)*OrienVecx(i,3,j)));
        
        Beta(i,1,j) = acosd((ConveyVec(1)*OrienVecy(i,1,j)+ ConveyVec(2)*OrienVecy(i,2,j) + ConveyVec(3)*OrienVecy(i,3,j)));
        
        Gamma(i,1,j) = acosd((ConveyVec(1)*OrienVecz(i,1,j)+ ConveyVec(2)*OrienVecz(i,2,j) + ConveyVec(3)*OrienVecz(i,3,j)));
        
    end
end

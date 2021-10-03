%% Orientation Quat Post Process
%-----------------------------------------------------

%Connor Murphy

%This code takes a dump file for multisphere that has quaternion
%and center of mass values:

% q0 q1 q2 q3 cmx cmy cmz


%-------------------------------------------------------
clear all
numPellets = 1
partsperpellet = 70
numTimesteps = 40000 ; 
dumpFrequency = 200 ;

DumpNumber = numTimesteps/dumpFrequency;

name = 'dump.pellet_orienrand';
fid = fopen(name,'r+');
Quat =  zeros(DumpNumber, 4, numPellets);
PosCM = zeros(DumpNumber, 3, numPellets);

for i = 1:DumpNumber
    for k = 1:9
        fgetl(fid);
    end
    for j = 1:numPellets
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

OrienVec = zeros(DumpNumber, 3, numPellets) 
x = 0; y = 0; z = 1;

for i = 1:DumpNumber
    for j = 1:numPellets
    q0 = Quat(i,1,j);
    q1 = Quat(i,2,j);
    q2 = Quat(i,3,j);
    q3 = Quat(i,4,j);
 
    Ans = [x*(q0^2+q1^2-q2^2-q3^2)+2*y*(q1*q2-q0*q3)+2*z*((q0*q2)+(q1*q3));
           2*x*(q0*q3+q1*q2) + y*(q0^2-q1^2+q2^2-q3^2)+ 2*z*((q2*q3)-(q0*q1));
           2*x*(q1*q3-q0*q2) + 2*y*(q0*q1 + q2*q3) + z*(q0^2-q1^2-q2^2+q3^2)]
    
    OrienVec(i,1,j)= Ans(1,1);
    OrienVec(i,2,j)= Ans(2,1);
    OrienVec(i,3,j)= Ans(3,1);
    end
end

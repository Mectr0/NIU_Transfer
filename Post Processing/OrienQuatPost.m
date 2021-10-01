%% Orientation Quat Post Process
%clear all
numPellets = 1
partsperpellet = 70
numTimesteps = 40000 ; 
dumpFrequency = 200 ;

DumpNumber = numTimesteps/dumpFrequency;

name = 'dump.pellet_orienrand';
fid = fopen(name,'r+');

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

for i = 1:DumpNumber
    q0 = Quat(i,1,1);
    q1 = Quat(i,2,1);
    q2 = Quat(i,3,1);
    q3 = Quat(i,4,1);
    %M = 2*[(q0^2+q1^2-0.5) (q1*q2-q0*q3) (q0*q2+q1*q3);
  %         (q0*q3+q1*q2) (q0^2+q2^2-0.5) (q2*q3-q0*q1);
  %         (q1*q3-q0*q2) (q0*q1+q2*q3) (q0^2+q3^2-0.5)];
    x = 0; y = 0; z = 1;
    %x = -.187; y = 0.8636; z = -0.4682; 
%    q0 = -0.137519;
%    q1 = -0.00769816;
%    q2 = 0.706887;
%    q3 = 0.693786;
    Ans = [x*(q0^2+q1^2-q2^2-q3^2)+2*y*(q1*q2-q0*q3)+2*z*((q0*q2)+(q1*q3));
           2*x*(q0*q3+q1*q2) + y*(q0^2-q1^2+q2^2-q3^2)+ 2*z*((q2*q3)-(q0*q1));
           2*x*(q1*q3-q0*q2) + 2*y*(q0*q1 + q2*q3) + z*(q0^2-q1^2-q2^2+q3^2)]
    %testing what unit vector i need to match orientation vector from xyz
    %dump
    AnsTest = [(q0^2+q1^2-q2^2-q3^2) 2*(q1*q2-q0*q3) 2*((q0*q2)+(q1*q3));
               2*(q0*q3+q1*q2) (q0^2-q1^2+q2^2-q3^2) 2*((q2*q3)-(q0*q1));
               2*(q1*q3-q0*q2)  2*(q0*q1 + q2*q3)  (q0^2-q1^2-q2^2+q3^2)];
   % F = [-.1870; 0.8636; -0.4682];
   F = [-.2053; .9787; 0]
    
    Test = AnsTest\F;
    
    OrienVec(i,1)= Ans(1,1);
    OrienVec(i,2)= Ans(2,1);
    OrienVec(i,3)= Ans(3,1);
end

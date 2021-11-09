%Orientation Post Process
clear all
close all
numPellets = 1;
n = 70 ;

numTimesteps = 40000 ; 
dumpFrequency = 200 ;

DumpNumber = numTimesteps/dumpFrequency;

TotalParticles = numPellets*n;
name = 'dump.posrand.xyz';
fid = fopen(name,'r+');
for i = 1:DumpNumber
   for k = 1:9
      fgetl(fid);
   end
   for j = 1:numPellets
     for ii = 1:n
         id = fscanf(fid,'%d',1);
         if ii == 1
           pos(i,1,ceil(id/n))= id;
           pos(i,2,ceil(id/n)) = fscanf(fid,'%g',1);
           pos(i,3,ceil(id/n)) = fscanf(fid,'%g',1);
           pos(i,4,ceil(id/n)) = fscanf(fid,'%g',1);
         else if ii == 8
           pos2(i,1,ceil(id/n)) = id;
           pos2(i,2,ceil(id/n)) = fscanf(fid,'%g',1);
           pos2(i,3,ceil(id/n)) = fscanf(fid,'%g',1);
           pos2(i,4,ceil(id/n)) = fscanf(fid,'%g',1);
         else
           fgetl(fid);
         end
         end
     end
   end
end

for i = 1:DumpNumber
    for j = 1:numPellets
      dx = pos2(i,2,j)-pos(i,2,j);
      dy = pos2(i,3,j)-pos(i,3,j);
      dz = pos2(i,4,j)-pos(i,4,j);
      Orien(i,1,j)= dx/norm([dx,dy,dz]);
      Orien(i,2,j)= dy/norm([dx,dy,dz]);
      Orien(i,3,j)= dz/norm([dx,dy,dz]);
      
    end
end
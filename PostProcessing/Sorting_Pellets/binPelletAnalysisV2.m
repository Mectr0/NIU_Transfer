% Function Created by Connor Murphy
% This function allows for addition of new pellets while a simulation runs.
function [pInEachBin, COM, Vel, Quat,nx,ny,nz] = binPelletAnalysisV2(name,nt,dx,dy,dz,maxP,n1)
fid = fopen(name,'r+');
for zz = 1:5
        fgetl(fid); % The dump files comes with 11 lines of script before data starts
end
%IMPORTANT: X and Z Bounds are hard-coded right now. These correspond to
%the size of the simulation domain
Xlo = 0.5;%fscanf(fid,'%g',1);
Xhi =   2;%fscanf(fid,'%g',1);
fscanf(fid,'%g',1);
fscanf(fid,'%g',1);
Ylo = fscanf(fid,'%g',1);
Yhi = fscanf(fid,'%g',1);
Zlo = 0;
fscanf(fid,'%g',1);
Zhi = 0.13;
fscanf(fid,'%g',1);
fgetl(fid);
fgetl(fid);


nx = ceil((Xhi - Xlo)/dx); %Number of bins in x direction
fprintf('Number of bins in X: %d  Bin size (m): %f \n', nx, dx)

ny = ceil((Yhi - Ylo)/dy); %Number of bins in y direction
fprintf('Number of bins in Y: %d  Bin size (m): %f \n', ny, dy)

nz = ceil((Zhi - Zlo)/dz); %Number of bins in z direction
fprintf('Number of bins in Z: %d  Bin size (m): %f \n', nz, dz)

%!!!!!!!!!NOTE I have to include domain outside of hopper unless I can
%disregard pellets that have left the hopper


pInEachBin = zeros(nx,ny,nz,nt,maxP); %Preallocation: maxP might need to be larger if bin size is increased

%Initial preallocation for different properties (3D for each timestep)
%IMPORTANT: First column for each value contain sorted pellet ids
COM = zeros(n1,4,nt); %COM: Center of mass [id, x, y, z]
Vel = zeros(n1,4,nt); %Vel: Velocity components [id, x, y, z]
Quat = zeros(n1,5,nt);%Quat: Quaternion values [id, q0, q1, q2, q3]

np = n1; % Initial number of pellets
for it = 1:nt %Loop through each timestep of the dump file
%----------------------------
    for zz = 1:3 % lines of script between each dumped timestep
        fgetl(fid);
    end
    
    npold = np; %Store a temporary copy of np
    np = fscanf(fid,'%d\n',1); %New number of pellets for current timestep
    if it == 1 && np ~= n1 %Making sure the initial number of pellets is correct
        error('n1 input must match initial number of pellets')
    end
    AddNum = np - npold; %The number of rows that must be added for newly inserted pellets
    
    counter=zeros(nx,ny,nz,nt); %Initialize the counter each timestep
    
    for zz = 1:5
        fgetl(fid);
    end
   
%-----------------------------
if AddNum > 0 %Concatenate new rows into the property matrices
    COM = [COM; zeros(AddNum,4,nt)];
    Vel = [Vel; zeros(AddNum,4,nt)];
    Quat =[Quat;zeros(AddNum,5,nt)]; 
end

   counter=zeros(nx,ny,nz,nt); %Reset the counter
   
   %Dump file: [id q0 q1 q2 q3 xx xy xz vx vy vz]
   tempcells = textscan(fid,'%f %f %f %f %f %f %f %f %f %f %f\n',np); %Read data from dump file
   tempmat = cell2mat(tempcells);
   
   %Sort data into proper matrices
   Quat(:,1:5,it) = tempmat(:,1:5); 
   COM(:,1:4,it) = tempmat(:,[1,6:8]);
   Vel(:,1:4,it) = tempmat(:,[1,9:11]);
   
   %Sort each matrix by pellet id      
   Quat(:,:,it) = sortrows(Quat(:,:,it)); 
   COM(:,:,it) = sortrows(COM(:,:,it));
   Vel(:,:,it) = sortrows(Vel(:,:,it));
   for ip = 1:np
      %Looping through pellets and finding which bin they belong to  
      i = ceil((COM(ip,2,it)-Xlo)/dx); %NOTE: DOMAIN MUST ALL BE POSITIVE
      j = ceil((COM(ip,3,it)-Ylo)/dy); 
      k = ceil((COM(ip,4,it)-Zlo)/dz); 
        
      %bID4p(ip,it,:) = [i,j,k];
      counter(i,j,k,it) = counter(i,j,k,it)+1; %Add a number to corresponding counter index to add pellets to the same bin
      pInEachBin(i,j,k,it,counter(i,j,k,it)) = COM(ip,1,it);
        
   end
end
fclose(fid);
end

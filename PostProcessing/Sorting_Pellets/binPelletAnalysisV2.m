% Function Created by Connor Murphy
function [pInEachBin, COM, Vel, Quat,nx,ny,nz] = binPelletAnalysisV2(name,nt,dx,dy,dz,maxP,n1)
fid = fopen(name,'r+');
for zz = 1:5
        fgetl(fid); % The dump files comes with 18 lines of script before data starts
end
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


nx = ceil((Xhi - Xlo)/dx); %ceil((Xhi - Xlo)/dx); NOTE: domain is large in this direction
fprintf('Number of bins in X: %d  Bin size (m): %f \n', nx, dx)

ny = ceil((Yhi - Ylo)/dy); %Number of bins in y direction
fprintf('Number of bins in Y: %d  Bin size (m): %f \n', ny, dy)

nz = ceil((Zhi - Zlo)/dz);
fprintf('Number of bins in Z: %d  Bin size (m): %f \n', nz, dz)

%!!!!!!!!!NOTE I have to include domain outside of hopper unless I can
%disregard pellets that have left the hopper

%bID4p = zeros(nmax,nt,3); %3D index
pInEachBin = zeros(nx,ny,nz,nt,maxP);


COM = zeros(n1,4,nt);
Vel = zeros(n1,4,nt);
Quat = zeros(n1,5,nt);
np = n1;
for it = 1:nt
%----------------------------
    for zz = 1:3 % lines of script between each dumped timestep
        fgetl(fid);
    end
    
    npold = np;
    np = fscanf(fid,'%d\n',1); %Number of pellets for this timestep
    if it == 1 && np ~= n1 %Making sure the initial number of pellets is correct
        error('n1 input must match initial number of pellets')
    end
    AddNum = np - npold;
    counter=zeros(nx,ny,nz,nt); %Initialize the counter each timestep
    
    for zz = 1:5
        fgetl(fid);
    end
   
%-----------------------------
if AddNum > 0
    COM = [COM; zeros(AddNum,4,nt)];
    Vel = [Vel; zeros(AddNum,4,nt)];
    Quat =[Quat;zeros(AddNum,5,nt)]; 
end

   counter=zeros(nx,ny,nz,nt); %Reset the counter
    
   tempcells = textscan(fid,'%f %f %f %f %f %f %f %f %f %f %f\n',np);
   tempmat = cell2mat(tempcells);
   Quat(:,1:5,it) = tempmat(:,1:5);
   COM(:,1:4,it) = tempmat(:,[1,6:8]);
   Vel(:,1:4,it) = tempmat(:,[1,9:11]);
        
        Quat(:,:,it) = sortrows(Quat(:,:,it)); %Sort each matrix by pellet id 
        COM(:,:,it) = sortrows(COM(:,:,it));
        Vel(:,:,it) = sortrows(Vel(:,:,it));
    for ip = 1:np
        
        i = ceil((COM(ip,2,it)-Xlo)/dx); %NOTE: DOMAIN MUST ALL BE POSITIVE
        j = ceil((COM(ip,3,it)-Ylo)/dy); 
        k = ceil((COM(ip,4,it)-Zlo)/dz); 
        
        %bID4p(ip,it,:) = [i,j,k];
        counter(i,j,k,it) = counter(i,j,k,it)+1;
        pInEachBin(i,j,k,it,counter(i,j,k,it)) = COM(ip,1,it);
        
    end
end
fclose(fid);
end

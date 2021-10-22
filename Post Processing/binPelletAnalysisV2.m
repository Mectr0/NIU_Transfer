function [pInEachBin, COM] = binPelletAnalysisV2(name,nt,nmax,dx,dy,dz,maxP)
fid = fopen(name,'r+');
for zz = 1:5
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


nx = ceil((2 - (0.5))/dx); %ceil((Xhi - Xlo)/dx); NOTE: domain is large in this direction
ny = ceil((Yhi - Ylo)/dy); %Number of bins in y direction
nz = ceil((Zhi - Zlo)/dz);

%!!!!!!!!!NOTE I have to include domain outside of hopper unless I can
%disregard pellets that have left the hopper

%bID4p = zeros(nmax,nt,3); %3D index
pInEachBin = zeros(nx,ny,nz,nt,maxP);
counter=zeros(nx,ny,nz,nt,maxP);

COM = zeros(nmax,4,nt);
Vel = zeros(nmax,4,nt);
Quat = zeros(nmax,5,nt);
for it = 1:nt
%----------------------------
    for zz = 1:3 % lines of script between each dumped timestep
        fgetl(fid);
    end
    np = fscanf(fid,'%d\n',1); %Number of pellets for this timestep
    
    for zz = 1:5
        fgetl(fid);
    end
   
%-----------------------------
    for ip = 1:np
        id = fscanf(fid,'%d',1);
        Quat(ip,1,it) = id;
        Quat(ip,2:5,it) = fscanf(fid,'%g %g %g %g',4);
        COM(ip,1,it) = id;
        COM(ip,2:4,it) = fscanf(fid,'%g %g %g',3);
        Vel(ip,1,it) = id;
        Vel(ip,2:4,it) = fscanf(fid,'%g %g %g\n',3);
        
        %!!!!!NOTE: I preallocate for max number of pellets. I need to have
        %a stop condition for empty rows in COM
    end
        Quat(:,:,it) = sortrows(Quat(:,:,it)); %Sort each matrix by pellet id 
        COM(:,:,it) = sortrows(COM(:,:,it));
        Vel(:,:,it) = sortrows(Vel(:,:,it));
    for ip = 1:np
        i = ceil(COM(ip,2)/dx); %NOTE: DOMAIN MUST ALL BE POSITIVE
        j = ceil(COM(ip,3)/dy);
        k = ceil(COM(ip,4)/dz); 
        %bID4p(ip,it,:) = [i,j,k];
        counter(i,j,k,it) = counter(i,j,k,it)+1;
        pInEachBin(i,j,k,it,counter(i,j,k,it)) = ip;
    end

end
fclose(fid)
end
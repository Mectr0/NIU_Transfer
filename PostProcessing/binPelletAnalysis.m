function [bID4p,pInEachBin] = binPelletAnalysis(name,nt,nmax,dx,dy,dz)
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
%!!!!!!!!!NOTE I have to include domain outside of hopper unless I can
%disregard pellets that have left the hopper
ny = ceil((Yhi - Ylo)/dy); %Number of bins in y direction
nz = ceil((Zhi - Zlo)/dz);

bID4p = zeros(nmax,nt,3); %3D index
pInEachBin = cell(nx,ny,nz,nt);
for it = 1:nt
%----------------------------
    for zz = 1:3 % lines of script between each dumped timestep
        fgetl(fid);
    end
    np = fscanf(fid,'%d\n',1); %Number of pellets for this timestep
    
    for zz = 1:5
        fgetl(fid);
    end
    Data = zeros(np,8);
%-----------------------------
    for ip = 1:np
        Data(ip,1:11) = fscanf(fid,'%d %g %g %g %g %g %g %g %g %g %g\n',11);%Pellet ID 
        %NOTE id q0 q1 q2 q3 xcm ycm zcm vxcm vycm vzcm
    end
    for ip = 1:np
        i = ceil(Data(ip,6)/dx); %NOTE DOMAIN MUST ALL BE POSITIVE
        j = ceil(Data(ip,7)/dy);
        k = ceil(Data(ip,8)/dz); 
        bID4p(ip,it,:) = [i,j,k];
        tmp = cell2mat(pInEachBin(i,j,k,it));
        pInEachBin(i,j,k,it) = mat2cell([tmp,Data(ip,1:11)],1);
    end

end
end
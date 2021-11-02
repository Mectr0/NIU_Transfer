function [pInEachBin, COM, Vel] = binPelletAnalysisV2(name,nt,dx,dy,dz,maxP,n1,px,py,pz)
fid = fopen(name,'r+');
for zz = 1:5
        fgetl(fid); % The dump files comes with 18 lines of script before data starts
end
Xlo = 0.4;%fscanf(fid,'%g',1);
Xhi = 2;%fscanf(fid,'%g',1);
fscanf(fid,'%g',1);
fscanf(fid,'%g',1);
Ylo = fscanf(fid,'%g',1);
Yhi = fscanf(fid,'%g',1);
Zlo = fscanf(fid,'%g',1);
Zhi = fscanf(fid,'%g',1);
fgetl(fid);
fgetl(fid);


nx = ceil((Xhi - Xlo)/dx); %ceil((Xhi - Xlo)/dx); NOTE: domain is large in this direction
ny = ceil((Yhi - Ylo)/dy); %Number of bins in y direction
nz = ceil((Zhi - Zlo)/dz);

%!!!!!!!!!NOTE I have to include domain outside of hopper unless I can
%disregard pellets that have left the hopper

%bID4p = zeros(nmax,nt,3); %3D index
pInEachBin = zeros(nx,ny,nz,nt,maxP);
counter=zeros(nx,ny,nz,nt,maxP);

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
    if it == 1 & np ~= n1 %Making sure the initial number of pellets is correct
        error('n1 input must match initial number of pellets')
    end
    AddNum = np - npold;
    
    for zz = 1:5
        fgetl(fid);
    end
   
%-----------------------------
    COM = [COM; zeros(AddNum,4,nt)];
    Vel = [Vel; zeros(AddNum,4,nt)];
    Quat =[Quat;zeros(AddNum,5,nt)]; 

    for ip = 1:np
        id = fscanf(fid,'%d',1);
        Quat(ip,1,it) = id;
        Quat(ip,2:5,it) = fscanf(fid,'%g %g %g %g',4);
        COM(ip,1,it) = id;
        COM(ip,2:4,it) = fscanf(fid,'%g %g %g',3);
        Vel(ip,1,it) = id;
        Vel(ip,2:4,it) = fscanf(fid,'%g %g %g\n',3);
    end
        Quat(:,:,it) = sortrows(Quat(:,:,it)); %Sort each matrix by pellet id 
        COM(:,:,it) = sortrows(COM(:,:,it));
        Vel(:,:,it) = sortrows(Vel(:,:,it));
    for ip = 1:np
%         if COM(ip,1) == 0 %Disregards empty rows
%         else
        i = ceil((COM(ip,2,it)-Xlo)/dx); %NOTE: DOMAIN MUST ALL BE POSITIVE
        j = ceil((COM(ip,3,it)-Ylo)/dy); %DOES THE DOMAIN HAVE TO START AT THE ORIGIN?
        k = ceil((COM(ip,4,it)-Zlo)/dz); 
        %if any(px(:) == i) && any(py(:) == j) && any(pz(:) == k)
        %bID4p(ip,it,:) = [i,j,k];
        counter(i,j,k,it) = counter(i,j,k,it)+1;
        pInEachBin(i,j,k,it,counter(i,j,k,it)) = COM(ip,1);
        %end
        end
    end
fclose(fid);
end

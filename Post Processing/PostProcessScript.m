clear
nmax = 300; %Max number of pellets
dumpfreq = 1000;
nsteps = 15000;
nt = nsteps/dumpfreq; %Number of timesteps
delta = 0.05; %0.006562; %Pellet Diameter : Should I use actual value or simulation value?
dx = 4*delta;
dy = 2.5*delta;
dz = 2*delta;
name = 'dump.pellet_orien10'; %Name of data file
maxP = 15 ;%Max amount of pellets that could fit in a bin (for preallocation)
n1 = 10; %number of pellets in the first timestep
px = [6:8];
py = [1:2];
pz = [1:2];

[pInEachBin,COM,Vel]= binPelletAnalysisV2(name,nt,dx,dy,dz,maxP,n1,px,py,pz);

for i = 1:nt

Inspect = nonzeros(pInEachBin(:,1:2,:,i,:));
VelavgY(i) = mean(Vel(Inspect,3,i));

end
plot(1:nt,VelavgY,'-ro')
xlabel('Timestep')
ylabel('average Y-velocity (m/s)')

%(nx,ny,nz,nt,maxP)

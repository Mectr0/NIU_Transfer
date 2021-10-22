clear
nmax = 300; %Max number of pellets
dumpfreq = 1000;
nsteps = 15000;
nt = nsteps/dumpfreq; %Number of timesteps
delta = 0.05; %0.006562; %Pellet Diameter : Should I use actual value or simulation value?
dx = 2*delta;
dy = 0.8*delta;
dz = 2*delta;
name = 'dump.pellet_orien'; %Name of data file
maxP = 15 ;%Max amount of pellets that could fit in a bin (for preallocation)

[pInEachBin,COM]= binPelletAnalysisV2(name,nt,nmax,dx,dy,dz,maxP);

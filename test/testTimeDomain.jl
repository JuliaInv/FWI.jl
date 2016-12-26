using  jInv.Mesh
using  jInv.Utils
using  jInv.LinearSolvers
using  jInv.InverseSolve
using  jInv.Vis
using  EikonalInv
using  MAT
using  FWI
using  ForwardHelmholtz
using  Multigrid
println("=======================================================================")
println("===============  TestTimeDomain   =====================================");
println("=======================================================================")

#############################################################################################################
modelDir = pwd();
dataDir = pwd();
include("../drivers/prepareFWIDataFiles.jl");
include("../drivers/readModelAndGenerateMeshMref.jl");

plotting = true;

if plotting
	using  PyPlot
	close("all");
end
#############################################################################################################
dataDir = pwd();
resultsDir = pwd();
########################################################################################################
dim     = 2;
pad     = 40;
ABLpad = pad+10;
jumpSrc = 80;
newSize = [250,120];
offset  = ceil(Int64,(newSize[1]*(13.5/13.5)));
println("Offset is: ",offset)
(m,Minv,mref,boundsHigh,boundsLow) = readModelAndGenerateMeshMref(modelDir,"SEGmodel2Dsalt.dat",dim,pad,[0.0,13.5,0.0,4.2],newSize,1.752,2.9);

omega = [1.5,2.0]*2*pi;

# ###################################################################################################################
dataFilenamePrefix = string(dataDir,"/DATA_SEG",tuple((Minv.n+1)...));
timeDataFilenamePrefix = string(dataDir,"/timeDATA_SEG",tuple((Minv.n+1)...));
resultsFilenamePrefix = string(resultsDir,"/SEG");
#######################################################################################################################

if plotting
	limits = [1.5,4.8];
	figure()
	plotModel(m,true,Minv,pad*0,limits);
end

# println("omega*maximum(h): ",omega*maximum(Minv.h)*sqrt(maximum(1./(boundsLow.^2))));


############################ Time Domain Data ##################################################
dt = 0.002;
T = 14.0;
fm = 3.0;
workerList = workers();
# prepareFWIDataFilesFromTime(m,Minv,mref,boundsHigh,boundsLow,timeDataFilenamePrefix,omega,pad,ABLpad,jumpSrc,offset,workerList,T,dt,fm,0.0000);

RCVfile = string(timeDataFilenamePrefix,"_rcvMap.dat");
SRCfile = string(timeDataFilenamePrefix,"_srcMap.dat");
srcNodeMap = readSrcRcvLocationFile(SRCfile,Minv);
rcvNodeMap = readSrcRcvLocationFile(RCVfile,Minv);
							
DobsTD = Array(Array{Complex128,2},length(omega));
WdTD = Array(Array{Complex128,2},length(omega));

for k = 1:length(omega)
	omRound = string(round((omega[k]/(2*pi))*100.0)/100.0);
	(Dk,Wk) =  readDataFileToDataMat(string(timeDataFilenamePrefix,"_freq",omRound,".dat"),srcNodeMap,rcvNodeMap);
	DobsTD[k] = Dk;
	WdTD[k] = Wk;
end
######################## Frequency Domain #################################################

file = matopen(string(timeDataFilenamePrefix,"_PARAM.mat"));
waveCoef = read(file,"waveCoef");
close(file);

Ainv = getJuliaSolver();

# This is a list of workers for FWI. Ideally they should be on different machines.
workersFWI = [workers()[1]];
println("The workers that we allocate for frequency FWI are:");
println(workersFWI)
maxBatchSize = 10;
# prepareFWIDataFiles(m,Minv,mref,boundsHigh,boundsLow,dataFilenamePrefix,omega,waveCoef, pad,ABLpad,jumpSrc,offset,workersFWI,maxBatchSize,Ainv);
### Read receivers and sources files
RCVfile = string(dataFilenamePrefix,"_rcvMap.dat");
SRCfile = string(dataFilenamePrefix,"_srcMap.dat");
srcNodeMap = readSrcRcvLocationFile(SRCfile,Minv);
rcvNodeMap = readSrcRcvLocationFile(RCVfile,Minv);

DobsFD = Array(Array{Complex128,2},length(omega));
WdFD = Array(Array{Complex128,2},length(omega));

for k = 1:length(omega)
	omRound = string(round((omega[k]/(2*pi))*100.0)/100.0);
	(Dk,Wk) =  readDataFileToDataMat(string(dataFilenamePrefix,"_freq",omRound,".dat"),srcNodeMap,rcvNodeMap);
	DobsFD[k] = Dk;
	WdFD[k] = Wk;
end   


# if plotting
	# for k=1:length(DobsTime)
		# figure();
		# imshow(DobsTime[k]);
	# end
# end



#######################################################
### Get Frequency domain data #########################
#######################################################
if plotting
	for jj = 1:length(omega) 
		figure(jj)
		subplot(1,3,1)
		imshow(real(DobsFD[jj])); colorbar()
		subplot(1,3,2)
		imshow(real(DobsTD[jj])); colorbar()
		subplot(1,3,3)
		imshow(abs(DobsTD[jj] - DobsFD[jj])); colorbar()
	end
end

##############################################################################################
# rm("DATA_SEG(120,60)_freq0.5.dat");
# rm("DATA_SEG(120,60)_freq1.0.dat");
# rm("DATA_SEG(120,60)_rcvMap.dat");
# rm("DATA_SEG(120,60)_srcMap.dat");
# rm("DATA_SEG(120,60)_PARAM.mat");
# rm("SEG(120,60)_FC2_HisGN.mat");
# rm("SEG(120,60)_FC1_HisGN.mat");
# rm("SEG(120,60)_Cyc1_FC2_HisGN.mat");
# rm("SEG(120,60)_Cyc1_FC1_HisGN.mat");
# rm("jInv.out");










   
   






   
   
   
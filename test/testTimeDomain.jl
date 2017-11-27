using  jInv.Mesh
using  jInv.Utils
using  jInv.LinearSolvers
using  jInv.InverseSolve
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

plotting = false;

realTesting = false;

if plotting
	using jInvVis
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

if realTesting
	jumpSrc = 35;
	newSize = [300,150];
else
	jumpSrc = 60;
	newSize = [200,100];
end
offset  = ceil(Int64,(newSize[1]*(8/13.5)));
println("Offset is: ",offset)
(m,Minv,mref,boundsHigh,boundsLow) = readModelAndGenerateMeshMref(modelDir,"SEGmodel2Dsalt.dat",dim,pad,[0.0,13.5,0.0,4.2],newSize,1.752,2.9);

if realTesting
	omega = [2.0,3.0]*2*pi;
else
	omega = [1.5,2.0]*2*pi;
end
	
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
if realTesting
	dt = 0.002;
	T  = 16.0;
	fm = 8.0;
else
	dt = 0.004;
	T = 10.0;
	fm = 3.0;
end


if plotting
	ricker,zeroOffset,widthOfHalfRicker = getRickerFunction(T,dt,fm);
	figure()
	plot(ricker);
end

workerList = workers();
DobsTimeNoisy = prepareFWIDataFilesFromTime(m,Minv,mref,boundsHigh,boundsLow,timeDataFilenamePrefix,omega,pad,ABLpad,jumpSrc,offset,workerList,T,dt,fm,0.005);

RCVfile = string(timeDataFilenamePrefix,"_rcvMap.dat");
SRCfile = string(timeDataFilenamePrefix,"_srcMap.dat");
srcNodeMap = readSrcRcvLocationFile(SRCfile,Minv);
rcvNodeMap = readSrcRcvLocationFile(RCVfile,Minv);
							
DobsTD = Array{Array{Complex128,2}}(length(omega));
WdTD = Array{Array{Complex128,2}}(length(omega));

for k = 1:length(omega)
	omRound = string(round((omega[k]/(2*pi))*100.0)/100.0);
	(Dk,Wk) =  readDataFileToDataMat(string(timeDataFilenamePrefix,"_freq",omRound,".dat"),srcNodeMap,rcvNodeMap);
	DobsTD[k] = Dk;
	WdTD[k] = Wk;
end


if plotting
	file = matopen(string(timeDataFilenamePrefix,"_timeDomain.mat"));
	DobsTime = read(file,"DobsTime");
	close(file);

	# for k=1:length(DobsTime)
		# figure();
		# imshow(DobsTime[k]);
	# end
end




DobsTTpicked = readDataFileToDataMat(string(timeDataFilenamePrefix,"_travelTime.dat"),srcNodeMap,rcvNodeMap)[1];

######################## Frequency Domain #################################################

# file = matopen(string(timeDataFilenamePrefix,"_PARAM.mat"));
# waveCoef = read(file,"waveCoef");
# close(file);

Ainv = getJuliaSolver();

# This is a list of workers for FWI. Ideally they should be on different machines.
workersFWI = [workers()[1]];
println("The workers that we allocate for frequency FWI are:");
println(workersFWI)
maxBatchSize = 10;
prepareFWIDataFiles(m,Minv,mref,boundsHigh,boundsLow,dataFilenamePrefix,omega,ones(Complex128,size(omega)), pad,ABLpad,jumpSrc,offset,workersFWI,maxBatchSize,Ainv,false,true);


### Read receivers and sources files
RCVfile = string(dataFilenamePrefix,"_rcvMap.dat");
SRCfile = string(dataFilenamePrefix,"_srcMap.dat");
srcNodeMap = readSrcRcvLocationFile(SRCfile,Minv);
rcvNodeMap = readSrcRcvLocationFile(RCVfile,Minv);

DobsFD = Array{Array{Complex128,2}}(length(omega));
WdFD = Array{Array{Complex128,2}}(length(omega));

for k = 1:length(omega)
	omRound = string(round((omega[k]/(2*pi))*100.0)/100.0);
	(Dk,Wk) =  readDataFileToDataMat(string(dataFilenamePrefix,"_freq",omRound,".dat"),srcNodeMap,rcvNodeMap);
	DobsFD[k] = Dk;
	WdFD[k] = Wk;
end   


for jj = 1:length(omega)
	println(vecnorm(DobsTD[jj] - DobsFD[jj]));
end


DobsTT = readDataFileToDataMat(string(dataFilenamePrefix,"_travelTime.dat"),srcNodeMap,rcvNodeMap)[1];



if plotting
	for jj = 1:length(omega) 
		figure()
		subplot(1,3,1)
		imshow(real(DobsFD[jj])); colorbar()
		subplot(1,3,2)
		imshow(real(DobsTD[jj])); colorbar()
		subplot(1,3,3)
		imshow(abs(DobsTD[jj] - DobsFD[jj])); colorbar()
	end
	figure()
	subplot(1,3,1)
	imshow(DobsTT); colorbar()
	subplot(1,3,2)
	imshow(DobsTTpicked); colorbar()
	subplot(1,3,3)
	imshow(abs(DobsTTpicked - DobsTT)); colorbar()
	
	file = matopen(string(timeDataFilenamePrefix,"_PARAM.mat"));
	zeroTimeOffset  = read(file,"zeroTimeOffset");
	dataMask 		= read(file,"dataMask");
	close(file);
	Tfirst = maximum(DobsTT) + 1.0;
	ifirst = round(Int64,Tfirst./dt);
	
	for k=1:length(DobsTime)
		Dtimek      = DobsTimeNoisy[k];
		mark = 1.5*maximum(abs(Dtimek));
		DTT_k       = DobsTT[:,k];
		DTTpicked_k = DobsTTpicked[:,k]
		for ii = 1:length(DTT_k)
			if dataMask[ii,k] == 1
				jj = convert(Int64,round(DTT_k[ii]./dt ) + zeroTimeOffset);
				Dtimek[ii,(jj-1):(jj+1)] = mark;
				jj = convert(Int64,round(DTTpicked_k[ii]./dt )+ zeroTimeOffset);
				Dtimek[ii,(jj-1):(jj+1)] = -mark;
			end
		end
		figure();
		subplot(2,1,1);
		imshow(Dtimek[:,1:ifirst]);colorbar() #imshow(dataTraces[k],clim=[0.005,-0.005]); colorbar();
		subplot(2,1,2);
		plot(vec(Dtimek[50,1:ifirst])+1)
		plot(vec(Dtimek[end-50,1:ifirst])-1.0,"r")
		plot(vec(Dtimek[div(size(Dtimek,1),2),1:ifirst]),"k")
	
		plot(vec(DobsTime[k][50,1:ifirst] +1.0),":")
		plot(vec(DobsTime[k][end-50,1:ifirst])-1.0,":r")
		plot(vec(DobsTime[k][div(size(Dtimek,1),2),1:ifirst]),":k")
	end
end

##############################################################################################

# error("Not deleting files!!!");

for k = 1:length(omega)
	omRound = string(round((omega[k]/(2*pi))*100.0)/100.0);
	rm(string(dataFilenamePrefix,"_freq",omRound,".dat"));
	timeDataFilenamePrefix
	rm(string(timeDataFilenamePrefix,"_freq",omRound,".dat"));
end

rm(string(timeDataFilenamePrefix,"_rcvMap.dat"));
rm(string(timeDataFilenamePrefix,"_srcMap.dat"));
rm(string(dataFilenamePrefix,"_rcvMap.dat"));
rm(string(dataFilenamePrefix,"_srcMap.dat"));
rm(string(dataFilenamePrefix,"_PARAM.mat"));
rm(string(timeDataFilenamePrefix,"_PARAM.mat"));
rm(string(timeDataFilenamePrefix,"_timeDomain.mat"));

rm(string(timeDataFilenamePrefix,"_travelTime.dat"));
rm(string(dataFilenamePrefix,"_travelTime.dat"));











   
   






   
   
   
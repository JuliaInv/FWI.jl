####
# This file contains the setup for running the 2D experiments presented in the paper
# Eran Treister and Eldad Haber, Full waveform inversion guided by 
# travel time tomography, 2017.
#
# The file should be run with a few workers already allocated. For example:
# Julia -p4 2017-TH-SISC_JointFwiEikDriverSEG.jl
####

using  jInv.Mesh
using  jInv.LinearSolvers
using  jInv.InverseSolve
using  jInv.Utils
using  jInvVis
using  EikonalInv
using  FWI
using  ForwardHelmholtz
using  Multigrid
using  MAT

experiment = "Joint";
# experiment = "TravelTimeInit"
# experiment = "OnlyFWI"

plotting = true;

if plotting
	using  PyPlot
	close("all");
end
#############################################################################################################
FWIDriversPath = 0;
modelDir = pwd();
dataDir = 0;
resultsDir = 0;

@everywhere FWIDriversPath = "../drivers/";
include(string(FWIDriversPath,"readModelAndGenerateMeshMref.jl"));
include(string(FWIDriversPath,"prepareFWIDataFiles.jl"));
include(string(FWIDriversPath,"setupJointInversion.jl"));
@everywhere include(string(FWIDriversPath,"remoteChangePmis.jl"));
dataDir = pwd();
timeDataDir = pwd();
resultsDir = pwd();


########################################################################################################
dim     = 2;
pad     = 30;
jumpSrc = 50;
newSize = [600,300];

offset  = ceil(Int64,(newSize[1]*(8.0/13.5)));
println("Offset is: ",offset," cells.")
(m,Minv,mref,boundsHigh,boundsLow) = readModelAndGenerateMeshMref(modelDir,"SEGmodel2Dsalt.dat",dim,pad,[0.0,13.5,0.0,4.2],newSize,1.752,2.9);


# omega = [2.0,2.5,3.5,4.5,6.0]*2*pi;
omega = [2.0,2.5]*2*pi;
maxBatchSize = 128;
useFilesForFields = false;


########################################################################################################
######################################## for 3D ########################################################
#######################################################################################################

# dim     	 = 3;
# pad     	 = 10;
# newSize 	 = [145,145,70];
# omega   	 = [1.5,2.0]*2*pi;
# jumpSrc 	 = 16;
# maxBatchSize     = 27;

# # calculation for num_sources^2: size = jumpSrc*(num_src-1) + 1 + 2*extraABL
# offset  = ceil(Int64,(newSize[1]*(8.0/13.5)));
# println("Offset is: ",offset," cells.")
# domain = [0.0,13.5,0.0,13.5,0.0,4.2];
# (m,Minv,mref,boundsHigh,boundsLow) = readModelAndGenerateMeshMref(modelDir,"3Dseg256256128.mat",dim,pad,domain,newSize,1.752,2.9);


# ###################################################################################################################
dataFilenamePrefix = string(dataDir,"/DATA_SEG",tuple((Minv.n+1)...));
timeDataFilenamePrefix = string(timeDataDir,"/TimeDATA_SEG",tuple((Minv.n+1)...));
resultsFilename = string(resultsDir,"/",experiment,"-TimeSEG");
#######################################################################################################################


writedlm(string(resultsFilename,tuple((Minv.n+1)...),"_mtrue.dat"),convert(Array{Float16},m));
writedlm(string(resultsFilename,tuple((Minv.n+1)...),"_mref.dat"),convert(Array{Float16},mref));

######################## ITERATIVE SOLVER FOR FWI #############################################
iterativeSolver = false;
Ainv = 0;
if iterativeSolver==true
	levels      = 3;
	numCores 	= 24;
	BLAS.set_num_threads(numCores);
	maxIter     = 50;
	relativeTol = 1e-4;
	relaxType   = "SPAI";
	relaxParam  = 1.0;
	relaxPre 	= 2;
	relaxPost   = 2;
	cycleType   ='W';
	coarseSolveType = "MUMPS";
	MG = getMGparam(levels,numCores,maxIter,relativeTol,relaxType,relaxParam,relaxPre,relaxPost,cycleType,coarseSolveType,0.0,0.0);
	shift 		= 0.15;
	Hparam = HelmholtzParam(Minv,zeros(0),zeros(0),0.0,true,true);
	Ainv = getShiftedLaplacianMultigridSolver(Hparam, MG,shift);
else   ######################## DIRECT SOLVER #################################################
	numCores 	= 16;
	BLAS.set_num_threads(numCores);
	Ainv = getMUMPSsolver([],0,0,2);
	# Ainv = getJuliaSolver();
end

##########################################################################################

println("omega*maximum(h): ",omega*maximum(Minv.h)*sqrt(maximum(1./(boundsLow.^2))));

ABLpad = pad + 4;

if dim==2
	# Generating frequency data from time tomain simulation - code works only for 2D.
	dt = 0.001;
	T  = 18.0;
	fm = 8.0;
	pickTraveltime = true;
	prepareFWIDataFilesFromTime(m,Minv,mref,boundsHigh,boundsLow,timeDataFilenamePrefix,omega,pad,ABLpad,jumpSrc,offset,workers(),T,dt,fm,0.01,pickTraveltime);
end

workersFWI = [workers()[1]];
println("The workers that we allocate for FWI are:");
println(workersFWI)
calcTravelTime = true;
prepareFWIDataFiles(m,Minv,mref,boundsHigh,boundsLow,dataFilenamePrefix,omega,ones(Complex128,size(omega)), pad,ABLpad,jumpSrc,offset,workersFWI,maxBatchSize,Ainv,useFilesForFields,calcTravelTime);

########################################################################################################################
################### READING AND COMPARING THE DATA FOR PLOTTING - NOT NECESSARY FOR INVERSION #######################################
########################################################################################################################
## Data that is generated through time domain (in 2D)
# if dim==2
	# RCVfile = string(timeDataFilenamePrefix,"_rcvMap.dat");
	# SRCfile = string(timeDataFilenamePrefix,"_srcMap.dat");
	# srcNodeMap = readSrcRcvLocationFile(SRCfile,Minv);
	# rcvNodeMap = readSrcRcvLocationFile(RCVfile,Minv);
							
	# DobsTD = Array(Array{Complex128,2},length(omega));
	# WdTD = Array(Array{Complex128,2},length(omega));

	# for k = 1:length(omega)
		# omRound = string(round((omega[k]/(2*pi))*100.0)/100.0);
		# (Dk,Wk) =  readDataFileToDataMat(string(timeDataFilenamePrefix,"_freq",omRound,".dat"),srcNodeMap,rcvNodeMap);
		# DobsTD[k] = Dk;
		# WdTD[k] = Wk;
	# end
	# DobsTTpicked = readDataFileToDataMat(string(timeDataFilenamePrefix,"_travelTime.dat"),srcNodeMap,rcvNodeMap)[1];
# end

## Data that is generated through frequency domain simulation
### Read receivers and sources files
# RCVfile = string(dataFilenamePrefix,"_rcvMap.dat");
# SRCfile = string(dataFilenamePrefix,"_srcMap.dat");
# srcNodeMap = readSrcRcvLocationFile(SRCfile,Minv);
# rcvNodeMap = readSrcRcvLocationFile(RCVfile,Minv);

# DobsFD = Array(Array{Complex128,2},length(omega));
# WdFD = Array(Array{Complex128,2},length(omega));

# for k = 1:length(omega)
	# omRound = string(round((omega[k]/(2*pi))*100.0)/100.0);
	# (Dk,Wk) =  readDataFileToDataMat(string(dataFilenamePrefix,"_freq",omRound,".dat"),srcNodeMap,rcvNodeMap);
	# DobsFD[k] = Dk;
	# WdFD[k] = Wk;
# end   

# DobsTT = readDataFileToDataMat(string(dataFilenamePrefix,"_travelTime.dat"),srcNodeMap,rcvNodeMap)[1];

########################################################################################################################
########################################################################################################################
########################################################################################################################
if dim == 2
	dataFilenamePrefix = timeDataFilenamePrefix;
end
(Q,P,pMis,SourcesSubInd,contDiv,Iact,sback,mref,boundsHigh,boundsLow,resultsFilename) = 
   setupJointInversion(m,dataFilenamePrefix,resultsFilename,plotting,workersFWI,maxBatchSize,Ainv,false,SSDFun,1.0,useFilesForFields);

########################################################################################################
# Setting up the inversion for slowness instead of velocity:
########################################################################################################
function dump(mc,Dc,iter,pInv,PMis,resultsFilename)
	fullMc = slowSquaredToVelocity(reshape(Iact*pInv.modelfun(mc)[1] + sback,tuple((pInv.MInv.n+1)...)))[1];
	Temp = splitext(resultsFilename);
	if iter>0
		Temp = string(Temp[1],iter,Temp[2]);
	else
		Temp = resultsFilename;
	end
	if resultsFilename!=""
		writedlm(Temp,convert(Array{Float16},fullMc));
	end
	if plotting
		figure(888,figsize = (22,10));
		clf();
		plotModel(fullMc,false,Minv,pad,[1.5,4.5],splitdir(Temp)[2]);
	end
end

#####################################################################################################
# Setting up the inversion for velocity:
#####################################################################################################
mref 		= velocityToSlowSquared(mref)[1];
t    		= copy(boundsLow);
boundsLow 	= velocityToSlowSquared(boundsHigh)[1];
boundsHigh 	= velocityToSlowSquared(t)[1]; t = 0;
modfun 		= identityMod;

EikPmis = contDiv[1]:contDiv[2]-1;
pMis[EikPmis] = multWd(pMis[EikPmis],50*sqrt(1./50.0));


########################################################################################################
# Set up Inversion #################################################################################
########################################################################################################

GN = "projGN"

maxStep=0.05*maximum(boundsHigh);

regparams = [1.0,1.0,1.0,1e-6];
# regfun(m,mref,M) = wdiffusionRegNodal(m,mref,M,Iact=Iact,C=regparams);
regfun(m,mref,M) = wFourthOrderSmoothingNodal(m,mref,M,Iact=Iact,C = regparams);

alpha = 1e+4;
pcgTol = 1e-1;
maxit = 15;
cgit = 5; 


if dim==2
	HesPrec=getExactSolveRegularizationPreconditioner();
else
	HesPrec = getSSORCGFourthOrderRegularizationPreconditioner(regparams,Minv,Iact,1.0,1e-8,1000);	
end 


pInv = getInverseParam(Minv,modfun,regfun,alpha,mref[:],boundsLow,boundsHigh,
                         maxStep=maxStep,pcgMaxIter=cgit,pcgTol=pcgTol,
						 minUpdate=1e-3, maxIter = maxit,HesPrec=HesPrec);
mc = copy(mref[:]);

EikPmis = contDiv[1]:contDiv[2]-1;
pMis[EikPmis] = multWd(pMis[EikPmis],50.0);

if experiment == "Joint"
	mode = "Joint";
	start = 2;
elseif experiment == "TravelTimeInit"
	mode = "1stInit";
	start = 1;
elseif experiment == "OnlyFWI"
	mode = "1stInit"
	start = 2;
end
mc,Dc = freqCont(mc, pInv, pMis,contDiv[1:3], 2, resultsFilename,dump,mode,start,1,GN);

### This is how we recover in case of a shut down in the middle of the run:
# mc = readdlm("OnlyFWI-TimeSEG(660,330)_Cyc1_FC2_GN15.dat");
# mc = Iact'*velocityToSlowSquared(mc[:])[1];

pInv.mref = copy(mc[:])
pMis[EikPmis] = multWd(pMis[EikPmis],sqrt(1./50.0));

regparams = [1.0,1.0,1.0,1e-6];
regfun2(m,mref,M) = wdiffusionRegNodal(m,mref,M,Iact=Iact,C=regparams);
pInv.regularizer = regfun2;
pInv.alpha = 1e+3;
pInv.maxIter = 5;
pInv.pcgMaxIter = 10;
if dim==3
	HesPrec = getSSORCGRegularizationPreconditioner(1.0,1e-5,1000);
	freqContBatch = 2;
else
	freqContBatch = 4;
end
pInv.HesPrec = HesPrec;

if experiment == "Joint"
	mode = "Joint";
	start = 2;
else
	mode = "1stInit"
	start = 2
end




mc,Dc = freqCont(mc, pInv, pMis,contDiv, freqContBatch, resultsFilename,dump,mode,start,2,GN);
pInv.alpha = 1e+2;
pInv.mref = mc[:]
mc,Dc = freqCont(mc, pInv, pMis,contDiv, freqContBatch, resultsFilename,dump,mode,start+1,3,GN);
pInv.alpha = 1e+1;
pInv.mref = mc[:]
mc,Dc = freqCont(mc, pInv, pMis,contDiv, freqContBatch, resultsFilename,dump,mode,start+1,4,GN);

#############################################################################################

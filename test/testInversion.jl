if nworkers()<=1
	addprocs(1);
end

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
println("=======================================================================\n")
println("===============  TestInversion using Julia solver   ===================\n");
println("=======================================================================\n")

#############################################################################################################
modelDir = pwd();


dataDir = pwd();
include("../drivers/prepareFWIDataFiles.jl");
include("../drivers/readModelAndGenerateMeshMref.jl");
include("../drivers/setupFWI.jl");

plotting = false;

if plotting
	using  PyPlot
	close("all");
end
#############################################################################################################
dataDir = pwd();
resultsDir = pwd();
########################################################################################################
dim     = 2;
pad     = 10;
jumpSrc = 15;
newSize = [100,50];
offset  = ceil(Int64,(newSize[1]*(13.5/13.5)));
println("Offset is: ",offset)
(m,Minv,mref,boundsHigh,boundsLow) = readModelAndGenerateMeshMref(modelDir,"SEGmodel2Dsalt.dat",dim,pad,[0.0,13.5,0.0,4.2],newSize,1.752,2.7);
maxBatchSize = 5;
omega = [0.5,1.0]*2*pi;

# ###################################################################################################################
dataFilenamePrefix = string(dataDir,"/DATA_SEG",tuple((Minv.n+1)...));
resultsFilenamePrefix = string(resultsDir,"/SEG");
#######################################################################################################################


if plotting
	limits = [1.5,4.8];
	figure()
	plotModel(m,true,true,Minv,pad*0,limits);
	figure()
	plotModel(mref,true,true,Minv,pad*0,limits);
end

######################## DIRECT SOLVER #################################################

# Ainv = getMUMPSsolver([],0,0,2);
Ainv = getJuliaSolver();
# Ainv = getPARsolver([],0,0,6);

##########################################################################################

println("omega*maximum(h): ",omega*maximum(Minv.h)*sqrt(maximum(1./(boundsLow.^2))));

# This is a list of workers for FWI. Ideally they should be on different machines.
workersFWI = [workers()[1]];
println("The workers that we allocate for FWI are:");
println(workersFWI)

prepareFWIDataFiles(m,Minv,mref,boundsHigh,boundsLow,dataFilenamePrefix,omega,one(Complex128)*ones(size(omega)), pad,pad + 10,jumpSrc,
					offset,workersFWI,maxBatchSize,Ainv);

(Q,P,pMis,SourcesSubInd,contDiv,Iact,sback,mref,boundsHigh,boundsLow,resultsFilename) = 
   setupFWI(m,dataFilenamePrefix,resultsFilenamePrefix,plotting,workersFWI,maxBatchSize,Ainv,SSDFun);

########################################################################################################
# Setting up the inversion for slowness
########################################################################################################
function dump(mc,Dc,iter,pInv,PMis,resultsFilename)
	fullMc = slowSquaredToVelocity(reshape(Iact*pInv.modelfun(mc)[1] + sback,tuple((pInv.MInv.n+1)...)))[1];
	# Temp = splitext(resultsFilename);
	# if iter>0
		# Temp = string(Temp[1],iter,Temp[2]);
	# else
		# Temp = resultsFilename;
	# end
	# if resultsFilename!=""
		# writedlm(Temp,convert(Array{Float16},fullMc));
	# end
	# if plotting
		# close(888);
		# figure(888);
		# plotModel(fullMc,true,false,[],0,[1.5,4.8],splitdir(Temp)[2]);
	# end
end

mref 		= velocityToSlow(mref)[1];
t    		= copy(boundsLow);
boundsLow 	= velocityToSlow(boundsHigh)[1];
boundsHigh 	= velocityToSlow(t)[1]; t = 0;
modfun 		= slowToSlowSquared;

########################################################################################################
# Set up Inversion #################################################################################
########################################################################################################

maxStep=0.05*maximum(boundsHigh);

regparams = [1.0,1.0,1.0,1e-4];
regfun(m,mref,M) = wdiffusionRegNodal(m,mref,M,Iact=Iact,C=regparams);
cgit = 2; 
alpha = 1e-10;
pcgTol = 1e-1;
maxit = 2;

HesPrec=getExactSolveRegularizationPreconditioner();

pInv = getInverseParam(Minv,modfun,regfun,alpha,mref[:],boundsLow,boundsHigh,
                         maxStep=maxStep,pcgMaxIter=cgit,pcgTol=pcgTol,
						 minUpdate=1e-3, maxIter = maxit,HesPrec=HesPrec);
mc = copy(mref[:]);
mc,Dc = freqCont(mc, pInv, pMis,contDiv, 3, resultsFilename,dump,"Joint",1,1,"projGN");

println("======================================================================================\n")
println("===============  TestInversion using shifted Laplacian multigrid   ===================\n")
println("======================================================================================\n")

######################## ITERATIVE SOLVER #############################################
levels      = 2;
numCores 	= 4;
BLAS.set_num_threads(numCores);
maxIter     = 50;
relativeTol = 1e-5;
relaxType   = "SPAI";
relaxParam  = 1.0;
relaxPre 	= 2;
relaxPost   = 2;
cycleType   ='W';
coarseSolveType = "NoMUMPS";
MG = getMGparam(levels,numCores,maxIter,relativeTol,relaxType,relaxParam,relaxPre,relaxPost,cycleType,coarseSolveType,0.0,0.0);
shift = 0.1;
AinvMG = getShiftedLaplacianMultigridSolver(Minv, MG,shift);


(Q,P,pMis,SourcesSubInd,contDiv,Iact,sback,mref,boundsHigh,boundsLow,resultsFilename) = 
   setupFWI(m,dataFilenamePrefix,resultsFilenamePrefix,plotting,workersFWI,maxBatchSize,AinvMG,SSDFun);

mref 		= velocityToSlow(mref)[1];
t    		= copy(boundsLow);
boundsLow 	= velocityToSlow(boundsHigh)[1];
boundsHigh 	= velocityToSlow(t)[1]; t = 0;
modfun 		= slowToSlowSquared;

pInv.mref = mref[:];					 
mc = copy(mref[:]);

mc,Dc = freqCont(mc, pInv, pMis,contDiv, 3, resultsFilename,dump,"Joint",1,0,"projGN");


##############################################################################################
rm("DATA_SEG(120,60)_freq0.5.dat");
rm("DATA_SEG(120,60)_freq1.0.dat");
rm("DATA_SEG(120,60)_rcvMap.dat");
rm("DATA_SEG(120,60)_srcMap.dat");
rm("DATA_SEG(120,60)_PARAM.mat");
rm("SEG(120,60)_FC2_HisGN.mat");
rm("SEG(120,60)_FC1_HisGN.mat");
rm("SEG(120,60)_Cyc1_FC2_HisGN.mat");
rm("SEG(120,60)_Cyc1_FC1_HisGN.mat");
rm("jInv.out");



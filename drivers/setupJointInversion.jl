function setupJointInversion(m,filenamePrefix::String,resultsOutputFolderAndPrefix::String,plotting::Bool,
		workersFWI::Array{Int64,1}=workers(),maxBatchSize::Int64 = 48,
		Ainv::AbstractSolver = getMUMPSsolver([],0,0,2),ignoreEik::Bool=false, misfun::Function=SSDFun,beta::Float64=1.0,useFilesForFields::Bool = false)
		
file = matopen(string(filenamePrefix,"_PARAM.mat"));
n_cells = read(file,"n");
OmegaDomain = read(file,"domain");
Minv = getRegularMesh(OmegaDomain,n_cells);
gamma = read(file,"gamma");
omega = read(file,"omega");
waveCoef = read(file,"waveCoef");
if length(omega)==1
	# matlab saves a 1 variable array as scalar.
	omega = [omega];
	waveCoef = [waveCoef];
end
warn("Assuming wavelet coeficient is taken into account in the observed data. RHS in FWI is just a delta func.");
waveCoef = ones(Complex128,size(waveCoef))

if ignoreEik==false
	HO = read(file,"HO");
end

boundsLow = read(file,"boundsLow");
boundsHigh = read(file,"boundsHigh");
mref =  read(file,"mref");
close(file);

resultsFilename = string(resultsOutputFolderAndPrefix,tuple((Minv.n+1)...),".dat");


### Read receivers and sources files
RCVfile = string(filenamePrefix,"_rcvMap.dat");
SRCfile = string(filenamePrefix,"_srcMap.dat");

srcNodeMap = readSrcRcvLocationFile(SRCfile,Minv);
rcvNodeMap = readSrcRcvLocationFile(RCVfile,Minv);

Q = generateSrcRcvProjOperators(Minv.n+1,srcNodeMap);
Q = Q.*1/(norm(Minv.h)^2);
println("We have ",size(Q,2)," sources");
P = generateSrcRcvProjOperators(Minv.n+1,rcvNodeMap);


########################################################################################################
##### Set up remote workers ############################################################################
########################################################################################################

N = prod(Minv.n+1);

Iact = speye(N);
mback   = zeros(Float64,N);
## Setting the sea constant:
mask = zeros(N);
mask[abs(m[:] .- minimum(m)) .< 1e-2] = 1.0;
mask[gamma[:] .>= 0.95*maximum(gamma)] = 1.0;
# setup active cells
mback = vec(mref[:].*mask);
# mback = vec(m[:].*mask);
sback = velocityToSlowSquared(mback)[1];
sback[mask .== 0.0] = 0.0;
Iact = Iact[:,mask .== 0.0];

boundsLow = Iact'*boundsLow;
boundsHigh = Iact'*boundsHigh;
mref = Iact'*mref[:];



####################################################################################################################
####################################################################################################################

if ignoreEik==false
	println("Reading Eik data:");	
	(DobsEik,WdEik) =  readDataFileToDataMat(string(filenamePrefix,"_travelTime.dat"),srcNodeMap,rcvNodeMap);

	EikMPIWorkers = nworkers(); # this just set the maximal MPI workers. To activate parallelism, run addprocs()

	(pForEik,contDivEik,SourcesSubIndEik) = getEikonalInvParam(Minv,Q,P,HO,EikMPIWorkers,useFilesForFields);
	Wd 		= Array(Array{Float64},length(pForEik));
	dobs 	= Array(Array{Float64},length(pForEik));
	for i=1:length(pForEik)
		I_i = SourcesSubIndEik[i];
		Wd[i]   = beta*WdEik[:,I_i];
		dobs[i] = DobsEik[:,I_i];
	end
	WdEik = 0;
	# DobsEik = 0;

	pMisEikRFs = getMisfitParam(pForEik, Wd, dobs, misfun, Iact,sback);
end
println("Reading FWI data:");

batch = min(size(Q,2),maxBatchSize);
(pForFWI,contDivFWI,SourcesSubIndFWI) = getFWIparam(omega,waveCoef,vec(gamma),Q,P,Minv,Ainv,workersFWI,batch,useFilesForFields);

# write data to remote workers
Wd   = Array(Array{Complex128,2},length(pForFWI))
dobs = Array(Array{Complex128,2},length(pForFWI))
for k = 1:length(omega)
	omRound = string(round((omega[k]/(2*pi))*100.0)/100.0);
	(DobsFWIwk,WdFWIwk) =  readDataFileToDataMat(string(filenamePrefix,"_freq",omRound,".dat"),srcNodeMap,rcvNodeMap);
	# figure()
	# imshow(real(DobsFWIwk)'); colorbar();title("Real")
	# figure()
	# imshow(imag(DobsFWIwk)'); colorbar();title("Imag")
	
	# figure()
	# imshow(real(exp(-1im*omega[k]*DobsEik).*DobsFWIwk)'); colorbar();
	# figure()
	# imshow(imag(exp(-1im*omega[k]*DobsEik).*DobsFWIwk)'); colorbar();
	
	# figure()
	# imshow(real(log(DobsFWIwk))'); colorbar(); title("Real log");
	# figure()
	# imshow(imag(log(DobsFWIwk))'); colorbar(); title("image log");
	
	for i = contDivFWI[k]:contDivFWI[k+1]-1
		I_i = SourcesSubIndFWI[i]; # subset of sources for ith worker.
		Wd[i] 	= WdFWIwk[:,I_i];
		# Wd[i]   = exp(-1im*omega[k]*DobsEik[:,I_i])
		dobs[i] = DobsFWIwk[:,I_i];
		
	end
	DobsFWIwk = 0;
	WdFWIwk = 0;
end


pMisFWIRFs = getMisfitParam(pForFWI, Wd, dobs, misfun, Iact,sback);

########################################################################################################
##### Set up remote workers ############################################################################
########################################################################################################

if ignoreEik==false
	SourcesSubInd    = [SourcesSubIndEik;SourcesSubIndFWI]; SourcesSubIndEik = 0;SourcesSubIndFWI = 0;
	pMis 			 = [pMisEikRFs;pMisFWIRFs];pMisEikRFs = 0;pMisFWIRFs = 0;
	contDiv			 = [contDivEik[1];contDivFWI + contDivEik[2] - 1];contDivEik = 0;contDivFWI = 0;
else
	SourcesSubInd    = SourcesSubIndFWI;
	pMis 			 = pMisFWIRFs;
	contDiv			 = contDivFWI;
end


return Q,P,pMis,SourcesSubInd,contDiv,Iact,sback,mref,boundsHigh,boundsLow,resultsFilename;
end

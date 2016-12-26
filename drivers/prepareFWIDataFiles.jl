function prepareFWIDataFiles(m,Minv::RegularMesh,mref,boundsHigh,boundsLow,
							filenamePrefix::String,omega::Array{Float64,1},waveCoef::Array{Complex128,1},
							pad::Int64,ABLpad::Int64,jump::Int64,offset::Int64=prod(Minv.n+1),workerList = workers(), 
							maxBatchSize::Int64=48, Ainv::AbstractSolver = getMUMPSsolver([],0,0,2),useFilesForFields::Bool = false)

########################## m is in Velocity here. ###################################
RCVfile = string(filenamePrefix,"_rcvMap.dat");
SRCfile = string(filenamePrefix,"_srcMap.dat");
writeSrcRcvLocFile(SRCfile,Minv,ABLpad,jump);
writeSrcRcvLocFile(RCVfile,Minv,ABLpad,1);
	
dataFullFilenamePrefix = string(filenamePrefix,"_freq");
gamma = prepareFWIDataFiles2(m, Minv, filenamePrefix,dataFullFilenamePrefix,omega,waveCoef,pad,ABLpad,offset,workerList,maxBatchSize,Ainv);	

file = matopen(string(filenamePrefix,"_PARAM.mat"), "w");
write(file,"boundsLow",boundsLow);
write(file,"boundsHigh",boundsHigh);
write(file,"mref",mref);
write(file,"domain",Minv.domain);
write(file,"n",Minv.n);
write(file,"gamma",gamma);
write(file,"pad",pad);
write(file,"omega",omega);
write(file,"waveCoef",waveCoef);
close(file);			
end

function prepareFWIDataFilesFromTime(m,Minv::RegularMesh,mref,boundsHigh,boundsLow,
							filenamePrefix::String,omega::Array{Float64,1},
							pad::Int64,ABLpad::Int64,jump::Int64,offset::Int64,workerList::Array{Int64},
							T::Float64,dt::Float64,fm::Float64,sigmaNoise::Float64 = 0.0)

########################## m is in Velocity here. ###################################
RCVfile = string(filenamePrefix,"_rcvMap.dat");
SRCfile = string(filenamePrefix,"_srcMap.dat");
writeSrcRcvLocFile(SRCfile,Minv,ABLpad,jump);
writeSrcRcvLocFile(RCVfile,Minv,ABLpad,1);
	
dataFullFilenamePrefix = string(filenamePrefix,"_freq");
gamma,rickt,waveCoef,Mask = prepareFWIDataFilesFromTime(m,Minv,filenamePrefix,dataFullFilenamePrefix,omega,pad,ABLpad,offset,workerList,T,dt,fm,sigmaNoise);

	

file = matopen(string(filenamePrefix,"_PARAM.mat"), "w");
write(file,"boundsLow",boundsLow);
write(file,"boundsHigh",boundsHigh);
write(file,"mref",mref);
write(file,"domain",Minv.domain);
write(file,"n",Minv.n);
write(file,"gamma",gamma);
write(file,"pad",pad);
write(file,"omega",omega);
write(file,"T",T);
write(file,"dt",dt);
write(file,"fm",fm);
write(file,"sigmaNoise",sigmaNoise);
write(file,"waveCoef",waveCoef);
write(file,"dataMask",Mask);
close(file);			
end




function prepareJointTravelTimeAndFWIDataFiles(m,Minv::RegularMesh,mref,boundsHigh,boundsLow,
							filenamePrefix::String,omega::Array{Float64,1},waveCoef::Array{Complex128,1},
							pad::Int64,ABLpad::Int64,jump::Int64,offset::Int64=prod(Minv.n+1),workerList = workers(), 
							maxBatchSize::Int64=48, Ainv::AbstractSolver = getMUMPSsolver([],0,0,2),useFilesForFields::Bool = false)

########################## m is in Velocity here. ###################################
RCVfile = string(filenamePrefix,"_rcvMap.dat");
SRCfile = string(filenamePrefix,"_srcMap.dat");
writeSrcRcvLocFile(SRCfile,Minv,ABLpad,jump);
writeSrcRcvLocFile(RCVfile,Minv,ABLpad,1);
	
dataFullFilenamePrefix = string(filenamePrefix,"_freq");
gamma = prepareFWIDataFiles2(m, Minv, filenamePrefix,dataFullFilenamePrefix,omega,waveCoef,pad,ABLpad,offset,workerList,maxBatchSize,Ainv,useFilesForFields);	
dataFullFilenamePrefix = string(filenamePrefix,"_travelTime");
HO = false;
prepareTravelTimeDataFiles(m, Minv, filenamePrefix,dataFullFilenamePrefix,offset,HO,useFilesForFields);

file = matopen(string(filenamePrefix,"_PARAM.mat"), "w");
write(file,"boundsLow",boundsLow);
write(file,"boundsHigh",boundsHigh);
write(file,"mref",mref);
write(file,"domain",Minv.domain);
write(file,"n",Minv.n);
write(file,"gamma",gamma);
write(file,"pad",pad);
write(file,"omega",omega);
write(file,"waveCoef",waveCoef);
write(file,"HO",HO);
close(file);			
end


function prepareFWIDataFiles2(m, Minv::RegularMesh, filenamePrefix::String,dataFullFilenamePrefix::String, omega::Array{Float64,1}, 
								waveCoef::Array{Complex128,1}, pad::Int64,ABLpad::Int64,offset::Int64,workerList::Array{Int64,1},maxBatchSize::Int64,
								Ainv::AbstractSolver,useFilesForFields::Bool = false)
########################## m is in Velocity here. ###################################
println("maxOmega*maxKappaSq*h: should be below 0.6");
println(omega[end]*maximum(Minv.h)*(1./minimum(m)));


RCVfile = string(filenamePrefix,"_rcvMap.dat");
SRCfile = string(filenamePrefix,"_srcMap.dat");
srcNodeMap = readSrcRcvLocationFile(SRCfile,Minv);
rcvNodeMap = readSrcRcvLocationFile(RCVfile,Minv);

Q = generateSrcRcvProjOperators(Minv.n+1,srcNodeMap);
P = generateSrcRcvProjOperators(Minv.n+1,rcvNodeMap);
Q = Q.*1/(norm(Minv.h)^2);
println("We have ",size(Q,2)," sources");
# compute observed data

ABLamp = getMaximalFrequency(1./(minimum(m).^2),Minv);
gamma = getABL(Minv,true,ones(Int64,Minv.dim)*ABLpad,ABLamp);
attenuation = 0.01*ABLamp;
gamma += attenuation; # adding Attenuation.

println("~~~~~~~ Getting data FWI: ~~~~~~~");
	
# Solve forward problem (should be relaced in reading data from file)

batch = min(size(Q,2),maxBatchSize);
(pFor,contDiv,SourcesSubInd) = getFWIparam(omega,waveCoef,vec(gamma),Q,P,Minv,Ainv,workerList,batch,useFilesForFields);

(D,pFor) = getData(velocityToSlowSquared(m[:])[1],pFor,ones(length(pFor)),true);

nsrc = size(Q,2);
nrcv = size(P,2);

for k = 1:length(omega)
	I = contDiv[k] : contDiv[k+1] - 1 
	Dobsk = Array(Array{Complex128,2},length(I));
	for i = 1:length(I)
		Dobsk[i] = fetch(D[I[i]]);
	end
	Dobsk = arrangeRemoteCallDataIntoLocalData(Dobsk);
	# Dobsk += 0.005*mean(abs(Dobsk))*(randn(size(Dobsk,1),size(Dobsk,2)) + 1im*randn(size(Dobsk,1),size(Dobsk,2)));
	omRound = string(round((omega[k]/(2*pi))*100.0)/100.0);
	Wd_k = (1./(abs(real(Dobsk))+0.1*mean(abs(Dobsk)))) + 1im*(1./(abs(imag(Dobsk))+0.1*mean(abs(Dobsk))));
	# Wd_k = (1./(0.0*abs(real(Dobsk))+1.0*mean(abs(Dobsk)))) + 1im*(1./(0.0*abs(imag(Dobsk))+1.0*mean(abs(Dobsk))));
	Wd_k = limitDataToOffset(Wd_k,srcNodeMap,rcvNodeMap,offset);
	filename = string(dataFullFilenamePrefix,omRound,".dat");
	writeDataFile(filename,Dobsk,Wd_k,srcNodeMap,rcvNodeMap);
end
return gamma;
end



function prepareFWIDataFilesFromTime(m, Minv::RegularMesh, filenamePrefix::String,dataFullFilenamePrefix::String, omega::Array{Float64,1}, 
									pad::Int64,ABLpad::Int64,offset::Int64,workerList::Array{Int64,1},T::Float64,dt::Float64,fm::Float64,sigmaNoise)
########################## m is in Velocity here. ###################################
println("maxOmega*maxKappaSq*h: should be below 0.6");
println(omega[end]*maximum(Minv.h)*(1./minimum(m)));


RCVfile = string(filenamePrefix,"_rcvMap.dat");
SRCfile = string(filenamePrefix,"_srcMap.dat");
srcNodeMap = readSrcRcvLocationFile(SRCfile,Minv);
rcvNodeMap = readSrcRcvLocationFile(RCVfile,Minv);

Q = generateSrcRcvProjOperators(Minv.n+1,srcNodeMap);
P = generateSrcRcvProjOperators(Minv.n+1,rcvNodeMap);
Q = Q.*1/(norm(Minv.h)^2);
println("We have ",size(Q,2)," sources");
# compute observed data


nsrc = size(Q,2);
nrcv = size(P,2);


ABLamp = getMaximalFrequency(1./(minimum(m).^2),Minv);
gamma = getABL(Minv,true,ones(Int64,Minv.dim)*ABLpad,ABLamp);
attenuation = 0.01*ABLamp;
gamma += attenuation; # adding Attenuation.

println("~~~~~~~ Getting time domain FWI data: ~~~~~~~");


Mask = ones(nrcv,nsrc);
Mask = limitDataToOffset(Mask,srcNodeMap,rcvNodeMap,offset);



rickt = getRickerFunction(T,dt,fm);
computeTimeDomain = false;

if computeTimeDomain
	pFor,SourcesSubInd = getTimeDomainFWIParam(gamma[:],Q,P,Mask,Minv,rickt,dt,T,workerList);
	tic()
	(DRF,pFor) = getData(velocityToSlowSquared(m[:])[1],pFor,ones(length(pFor)),true);
	toc()
	DobsTime = Array(Array{Float64,2},size(Q,2));
	for k=1:length(pFor)
		Dt = fetch(DRF[k]);
		for i = 1:length(SourcesSubInd[k])
			DobsTime[SourcesSubInd[k][i]] = Dt[i];
		end
		Dt = 0;
	end
	file = matopen(string(filenamePrefix,"_timeDomain.mat"),"w");
	write(file,"DobsTime",convert(Array{Array{Float32,2},1},DobsTime));
	close(file);
else
	file = matopen(string(filenamePrefix,"_timeDomain.mat"));
	DobsTime = read(file,"DobsTime");
	close(file);
end


DobsOmHat = Array(Array{Complex128,2},length(omega));
for k=1:length(omega)
	DobsOmHat[k] = zeros(Complex128,nrcv,nsrc);
end

nw = length(rickt);
rickh = fft(rickt);
rickhsub = zeros(Complex128,length(omega));

for k = 1:length(DobsTime)
	Dk = DobsTime[k];
	Dk = Dk + sigmaNoise*maximum(abs(Dk))*randn(size(Dk));
	Dhat = fft(Dk,2);
	for jj = 1:length(omega)
		iw = (nw*dt*omega[jj])/(2*pi) + 1
		iw = convert(Int64,round(iw));
		DobsOmHat[jj][:,k] = Dhat[:,iw];
		rickhsub[jj] = rickh[iw];
	end
end
for k = 1:length(omega)
	omRound = string(round((omega[k]/(2*pi))*100.0)/100.0);
	Dobsk = DobsOmHat[k];
	Wd_k = (1./(abs(real(Dobsk))+0.1*mean(abs(Dobsk)))) + 1im*(1./(abs(imag(Dobsk))+0.1*mean(abs(Dobsk))));
	Wd_k .*= Mask;
	filename = string(dataFullFilenamePrefix,omRound,".dat");
	writeDataFile(filename,Dobsk,Wd_k,srcNodeMap,rcvNodeMap);
end
return gamma,rickt,rickhsub,Mask;
end









using FWI
using jInv.Mesh;
using jInv.LinearSolvers
using jInv.InverseSolve
using jInv.Utils
using EikonalInv
using ForwardHelmholtz
using Multigrid

include("getAnalyticalMediums.jl");

function checkDerivativeFWI(f::Function,x0,param,Iact;out::Bool=true,tol::Number=1.9,nSuccess::Int=3)
	if out
		println(@sprintf("%9s\t%9s\t%9s\t%9s\t%9s\t%5s","h","E0","E1","O1","O2","OK?"))
	end
	v = 0.1*randn(size(x0))
	v = Iact*Iact'*v;
	# if eltype(x0)<:Complex
		# v += 1im*randn(size(v))
	# end
	f0,dvf = f(x0,v,param)
	
	
	# dvf   = real(dvf)
	Error = zeros(10,2)
	Order = zeros(10,2)
	Success = zeros(10)
	for j=1:5
		ft = f(x0+10.0^(-j)*v,zeros(size(x0)),param)[1]            # function value
		Error[j,1] = norm(f0-ft)/norm(f0)           # Error TaylorPoly 0
		Error[j,2] = norm(f0 .+10.0^(-j)*dvf .- ft)/norm(f0) # Error TaylorPoly 1
		if j>1
			Order[j,:] = log10(Error[j-1,:]./Error[j,:]);
		end
		if (Order[j,2]>tol) || (Error[j,1]/Error[j,2] > 100); Success[j]=1; end
		if out 
			println(@sprintf("%1.3e\t%1.3e\t%1.3e\t%1.3e\t%1.3e\t%5d",
							10.0^(-j), Error[j,1],Error[j,2], Order[j,1],Order[j,2],Success[j]))
		end
	end
	pass = sum(Success) > nSuccess
	return  pass,Error,Order
end


println("########### 2D TEST SETTING ##############################");
n = [50,100]
Omega = [0.0,4.0,0.0,8.0];
Minv = getRegularMesh(Omega,n-1);
println(Minv.n)
m = getAnalyticalConstGradInv2D(Minv.n+1,Minv.h)[1];
n = [100,50]
Omega = [0.0,8.0,0.0,4.0];
Minv = getRegularMesh(Omega,n-1);
m = m';

# using PyPlot
# imshow(m')
pad = 16;
jump = 15;
offset = 100;

# omega = [0.5./maximum(Minv.h), 0.2./maximum(Minv.h)];
omega = [0.5./maximum(Minv.h)];

WaveletCoef = rand(Complex128,size(omega));
println(omega)


RCVfile = string("rcvMap.dat");
SRCfile = string("srcMap.dat");
writeSrcRcvLocFile(SRCfile,Minv,pad,jump);
writeSrcRcvLocFile(RCVfile,Minv,pad,2);
srcNodeMap = readSrcRcvLocationFile(SRCfile,Minv);
rcvNodeMap = readSrcRcvLocationFile(RCVfile,Minv);
Q = generateSrcRcvProjOperators(Minv.n+1,srcNodeMap);
P = generateSrcRcvProjOperators(Minv.n+1,rcvNodeMap);

ABLPad = pad+5;
gamma = getABL(Minv,true,ones(Int64,Minv.dim)*ABLPad,omega[1]);
attenuation = 0.01*omega[1];
gamma += attenuation; # adding Attenuation.

# Ainv = getMUMPSsolver([],0,0,0);
# Ainv = getJuliaSolver();

###################################################################################################
levels      = 2;
numCores 	= 4;
BLAS.set_num_threads(numCores);
maxIter     = 50;
relativeTol = 1e-4;
relaxType   = "SPAI";
relaxParam  = 1.0;
relaxPre 	= 2;
relaxPost   = 2;
cycleType   ='W';
coarseSolveType = "NoMUMPS";
MG = getMGparam(levels,numCores,maxIter,relativeTol,relaxType,relaxParam,relaxPre,relaxPost,cycleType,coarseSolveType,0.0,0.0);
shift = 0.15;
Ainv = getShiftedLaplacianMultigridSolver(Minv, MG,shift);


###################################################################################################



println("Number of Sources: ",size(Q,2));

(pFor,contDiv,SourcesSubInd) = getFWIparam(omega,WaveletCoef,vec(gamma),Q,P,Minv,Ainv,workers(),16);

# setup active cells - this is important so that Sommerfeld BC will not interfere with the sensitivities.

N = prod(Minv.n+1);
Iact = speye(N);
mask = zeros(N);
mask[gamma[:] .>= 0.2*maximum(gamma)] = 1;
Iact = Iact[:,mask .== 0.0];


function f(m,v=zeros(size(m)),pFor = [])
	(D,pFor) = getData(m[:],pFor);
	D = fetch(D[1]);
	pf = fetch(pFor[1]);
	Jv = 0;
	if vecnorm(v)>0.0
		Jv = getSensMatVec(v[:],m[:],pf);
	end
	return D[:],Jv;
end
checkDerivativeFWI(f,m[:],pFor,Iact);

pf = fetch(pFor[1]);
v = rand(N);
Jv = getSensMatVec(v,m[:],pf);
Z = rand(Complex128,size(Jv));
JtZ = getSensTMatVec(Z,m[:],pf);

println(dot(Jv[:],Z[:]))
println(dot(JtZ[:],v))

if real(dot(JtZ[:],v) - dot(Jv[:],Z[:])) > 1e-5
	warn("Checking that (J^T)^T = J for FWI falied: ",real(dot(JtZ[:],v) - dot(Jv[:],Z[:])) );
end

############################## STOCHASTIC #####################################################
# println("####################### STOCHASTIC SELECTION ################################")
# sourcesSelection = randperm(size(Q,2))[1:div(size(Q,2),2)];
# setSourceSelection(pFor[1], sourcesSelection)

# function f(m,v=zeros(size(m)),pFor = [])
	# (D,pFor) = getData(m[:],pFor,Mesh2Mesh);
	# D = fetch(D[1]);
	# pf = fetch(pFor[1]);
	# Jv = 0;
	# if vecnorm(v)>0.0
		# Jv = getSensMatVec(v[:],m[:],pf);
	# end
	# return D[:],Jv;
# end
# checkDerivativeFWI(f,m[:],pFor,Iact);

# pf = fetch(pFor[1]);
# v = rand(N);
# Jv = getSensMatVec(v,m[:],pf);
# Z = rand(Complex128,size(Jv));
# JtZ = getSensTMatVec(Z,m[:],pf);

# println(dot(Jv[:],Z[:]))
# println(dot(JtZ[:],v))

# if real(dot(JtZ[:],v) - dot(Jv[:],Z[:])) > 1e-8
	# warn("Checking that (J^T)^T = J for FWI falied: ",real(dot(JtZ[:],v) - dot(Jv[:],Z[:])) );
# end


rm("srcMap.dat");
rm("rcvMap.dat");



export timeDomainFWIparam, simulateTimeDomainFWI,getTimeDomainFWIParam,getData,getRickerFunction

type timeDomainFWIparam 
    gamma					:: Vector{Float64}     # attenuation
    Sources					:: Union{Vector{Float64},SparseMatrixCSC,Array{Float64,2}}   # Sources
    Receivers				:: Union{Vector{Float64},SparseMatrixCSC,Array{Float64,2}}
	ReceiverMask			:: Array
	Mesh      				:: RegularMesh
	ricker					:: Array
	dt						:: Float64
	T 						:: Float64
end

function getTimeDomainFWIParam(gamma,Sources,Receivers,ReceiverMask,Mesh,ricker,dt,T,workerList::Array{Int64}=[])
if workerList==[]
	ActualWorkers = workers();
else
	ActualWorkers = intersect(workerList,workers());
	if length(ActualWorkers)<length(workerList)
		warn("FWI: workerList included indices of non-existing workers.")
	end
end
numWorkers = length(ActualWorkers);
println("We have ",numWorkers," workers for time FWI");
pFor   = Array(RemoteChannel,numWorkers);
SourcesSubInd = Array(Array{Int64,1},numWorkers);
i = 1; nextidx() = (idx=i; i+=1; idx)
nsrc  = size(Sources,2);

# send out jobs
@sync begin
	for p=ActualWorkers
		@async begin
			while true
				idx = nextidx()
				if idx > numWorkers 
					break
				end
				I_k = getSourcesIndicesOfKthWorker(numWorkers,idx,nsrc)
				SourcesSubInd[idx] = I_k;
				pFor[idx] = initRemoteChannel(timeDomainFWIparam,p, gamma,Sources[:,I_k], Receivers, ReceiverMask[:,I_k], Mesh, ricker,dt,T);
				wait(pFor[idx]);
			end
		end
	end
end
return pFor, SourcesSubInd# Array of Remote Refs
end

function getRickerFunction(T,dt,fm)
	Td = 2*sqrt(6)/(pi*fm);
	nt = convert(Int64,round(T/dt));
	t  = dt:dt:nt*dt;
	ricker = (1-(2*pi^2*fm^2).*((t-3*Td).^2)).*exp(-((t-3*Td).^2*(pi*fm).^2));
	return ricker;
end


function getData(m,pFor::timeDomainFWIparam,doClear::Bool=false)
Mesh       	= pFor.Mesh
gamma   	= pFor.gamma
Q       	= pFor.Sources
P       	= pFor.Receivers
dt 			= pFor.dt
T           = pFor.T;
rickFunc    = pFor.ricker
RecMask     = pFor.ReceiverMask
	
nrec  		= size(P,2) 
nsrc  		= size(Q,2)
	
D = Array(Array{Float64,2},nsrc);
	
nt = convert(Int64,round(T/dt));
# A = ForwardHelmholtz.getNodalLaplacianMatrix(Mesh);

N = prod(Mesh.n+1);
U = zeros(Float64,N,nt);


temp = zeros(N);
## m is assumed to be slowness squared here.
vsq = 1./m;
gamma = gamma.*(dt/2);

for k = 1:nsrc
	# Dk = zeros(Float64,sum(RecMask[:,k] != 0),nt);
	q = vec(Q[:,k]).*vsq;
	# println("Simulating source ",k," in time domain");
	# tic()
	for i = 2:nt-1
		# U[:,i+1] = (-dt^2*(v.^2).*(A*U[:,i]) + 2*U[:,i] - U[:,i-1] + gamma.*(dt/2).*U[:,i-1]+ rickFunc[i]*dt^2*((v.^2).*q)  )./(1+gamma*(dt/2));
		multOpNeumann!(Mesh,U[:,i],temp,Lap2DStencil);
		# U[:,i+1] = (-dt^2*vsq.*temp + 2*U[:,i] - U[:,i-1] + gamma.*U[:,i-1]+ (rickFunc[i]*dt^2)*q)./(1+gamma);
		@simd for j=1:size(U,1)
			@inbounds U[j,i+1] = (-dt^2*vsq[j].*temp[j] + 2*U[j,i] - U[j,i-1] + gamma[j].*U[j,i-1]+ (rickFunc[i]*dt^2)*q[j])./(1+gamma[j]);
		end
	end
	# toc()
	Dk = P'*U;
	Dk = vec(RecMask[:,k].!=0.0).*Dk;
	D[k] = Dk;
	U[:] = 0.0;
end
return D,pFor
end
export timeDomainFWIparam, simulateTimeDomainFWI,getTimeDomainFWIParam,getData,getRickerFunction

type timeDomainFWIparam <: ForwardProbType
    gamma					:: Vector{Float64}     # attenuation
    Sources					:: Union{Vector{Float64},SparseMatrixCSC,Array{Float64,2}}   # Sources
    Receivers				:: Union{Vector{Float64},SparseMatrixCSC,Array{Float64,2}}
	ReceiverMask			:: Array{Int8}
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
pFor   = Array{RemoteChannel}(numWorkers);
SourcesSubInd = Array{Array{Int64,1}}(numWorkers);
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
	ricker = (1-(2*pi^2*fm^2).*((t-3*Td).^2)).*exp.(-((t-3*Td).^2*(pi*fm).^2));
	zeroOffset = convert(Int64,round(((3*Td)./dt) + 1));
	widthOfHalfRicker = round(Int64,(sqrt(8)/(pi*fm))/dt);
	return ricker,zeroOffset,widthOfHalfRicker;
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
	
D = Array{Array{Float32,2}}(nsrc);
	
nt = convert(Int64,round(T/dt));

# A = ForwardHelmholtz.getNodalLaplacianMatrix(Mesh);

N = prod(Mesh.n+1);

u_ip1 = zeros(N);
u_im1 = zeros(N);
u_i   = zeros(N);
temp = zeros(N);

## m is assumed to be slowness squared here.
vsq = dt^2./m;
gammadt = gamma.*(0.5*dt);
tic();
for k = 1:nsrc
	println("At source ",k," out of ",nsrc);
	q = vec(Q[:,k]).*vsq;
	Dk = calcTrace(Mesh,N,nt,P,vsq,gammadt,q,u_ip1,u_im1,u_i,temp,rickFunc);
	D[k] = vec(RecMask[:,k]).*Dk;
	u_ip1[:] = 0.0;
	u_im1[:] = 0.0;
	u_i[:] = 0.0;
end
println("time domain get data took ",toc()," seconds.");

return D,pFor
end

function calcTrace(Mesh::RegularMesh,n::Int64,nt::Int64,P::SparseMatrixCSC,vsq::Array{Float64,1},gammadt::Array{Float64,1},q::Array{Float64,1},u_ip1::Array{Float64,1},u_im1::Array{Float64,1},u_i::Array{Float64,1},temp::Array{Float64,1},rickFunc::Array{Float64})
Dk    = zeros(Float32,size(P,2),nt);
# tic()	
# for i = 2:nt-1
	## Code with 3 u's
	# multOpNeumann!(Mesh,u_i,temp,Lap2DStencil);
	# u_ip1 = (-vsq.*temp + 2.0*u_i - u_im1 + gammadt.*u_im1+ rickFunc[i]*q)./(1.0+gammadt);
	# Dk[:,i+1] = P'*u_ip1;
	# t = u_im1;
	# u_im1 = u_i;
	# u_i = u_ip1;
	# u_ip1 = t;
# end
# toc()

for i = 2:nt-1
	### Code with 3 u's
	multOpNeumann!(Mesh,u_i,temp,Lap2DStencil);
	rfi = rickFunc[i];
	@simd for j=1:n
		@inbounds temp[j] = -vsq[j]*temp[j] + rfi*q[j]
	end
	@simd for j=1:n
		@inbounds u_ip1[j] = (temp[j] + 2.0*u_i[j] - u_im1[j] + gammadt[j]*u_im1[j] )/(1.0+gammadt[j]);
	end
	Dk[:,i+1] = P'*u_ip1;
	t = u_im1;
	u_im1 = u_i;
	u_i = u_ip1;
	u_ip1 = t;
end

return Dk;
end
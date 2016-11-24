
if nworkers()<2
	addprocs(3-nworkers())
end

using jInv.Mesh
using jInv.LinearSolvers
using jInv.ForwardShare
using ForwardHelmholtz
using FWI
using Base.Test


# get a small mesh
domain = [0 2.0 0 2.2]
n      = [12;8;]
M      = getRegularMesh(domain,n)
omega  = 0.1; 

# sources / receivers on top edge
idx    = reshape(collect(1:prod(n+1)),tuple(n+1...))
ib     = idx[:,1];
I      = speye(prod(n+1))
Q      = I[:,vec(ib)]
R      = copy(Q)

# get param without parallelizaion
ABLpad = 2;
gamma = getABL(M,true,ones(Int64,M.dim)*ABLpad,1.0);
attenuation = 0.01;
gamma += attenuation; # adding Attenuation.
gamma = gamma[:];
Ainv = getJuliaSolver();

pFor = getFWIparam(omega, one(Complex128),gamma,Q,R,M,Ainv,[workers()[1]])[1];

m0   = rand(tuple(n+1...))+1.0
dobs, = getData(vec(m0),fetch(pFor[1]));


# parallelize over sources
pForp,continuationDivision,SourcesSubInd = getFWIparam(omega, one(Complex128),gamma,Q,R,M,Ainv,workers())
dobsRF, = getData(vec(m0),pForp)
Dobs = zeros(eltype(dobs),size(dobs));
for k=1:length(dobsRF)
	dobsk = fetch(dobsRF[k]);
	Dobs[:,SourcesSubInd[k]] = dobsk;
end

@test vecnorm(Dobs-dobs)/vecnorm(dobs) < 1e-8




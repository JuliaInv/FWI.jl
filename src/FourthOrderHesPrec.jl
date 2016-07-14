
export getSSORCGFourthOrderRegularizationPreconditioner,SSORCGParamFourthOrder, applySSORCGFourth, setupDoNothing, wFourthOrderSmoothingNodal


function wFourthOrderSmoothingNodal(m::Vector, mref::Vector, M::AbstractMesh; Iact=1.0, C=[])	
	dm = m.-mref;
	d2R = wdiffusionRegNodal(m,mref,M,Iact = Iact,C = C)[3];
	clear(M);
	d2R = d2R'*d2R;
	dR  = d2R*dm;
	Rc  = 0.5*dot(dm,dR);
   return Rc,dR,d2R
end



import jInv.InverseSolve.HessianPreconditioner
type SSORCGParamFourthOrder
	C 			::Vector;
	M			::AbstractMesh;
	Iact			
	omega		::Float64
	tol			::Float64
	maxCGIter	::Int64
end

function getSSORCGFourthOrderRegularizationPreconditioner(C,M::AbstractMesh,Iact, omega::Float64=1.0,tol::Float64=1e-2,maxCGIter::Int64=100)
	return HessianPreconditioner(SSORCGParamFourthOrder(C,M,Iact,omega,tol,maxCGIter),applySSORCGFourth,setupDoNothing);
end

function getSAAMGFourthOrderRegularizationPreconditioner(C,M::AbstractMesh,Iact, omega::Float64=1.0,tol::Float64=1e-2,maxCGIter::Int64=100)
	return HessianPreconditioner(SSORCGParamFourthOrder(C,M,Iact,omega,tol,maxCGIter),applySAAMGFourth,setupDoNothing);
end


function applySAAMGFourth(Hs::Function, d2R::SparseMatrixCSC,v::Vector,param)
	d2R = wdiffusionRegNodal(v,v,param.M,Iact = param.Iact, C = param.C)[3]
	MG = getMGparam(5,8,500,1e-8,"SPAI",1.0,1,1,'V',"MUMPS");
	SA_AMGsetup(d2R,MG,eltype(v),true,1,false);
	x = zeros(eltype(v),size(v));
	x = solveCG_MG(d2R,MG,v,x,false)[1];
	v = copy(x); x[:] = 0.0;
	MG.relativeTol = 1e-4;
	x = solveCG_MG(d2R,MG,v,x,false)[1];
	return x;
end

function applySSORCGFourth(Hs::Function, d2R::SparseMatrixCSC,v::Vector,param)
	d2R = wdiffusionRegNodal(v,v,param.M,Iact = param.Iact, C = param.C)[3]
	Dinv = param.omega./diag(d2R);
	aux   = zeros(size(d2R,2));
	SSOR(r) = (aux[:]=0.0; ssorPrecTrans!(d2R,aux,r,Dinv); return aux);
	# println("In prec.")
	# tic();
	x = KrylovMethods.cg(d2R,v,tol=param.tol,maxIter=param.maxCGIter,M=SSOR,out=-1)[1]
	x = KrylovMethods.cg(d2R',copy(x),tol=param.tol,maxIter=param.maxCGIter,M=SSOR,out=-1)[1]
	# println("Done prec in ",toq()," seconds.");
	return x;
end

function setupDoNothing(Hs::Function, d2R::SparseMatrixCSC,param)
	return;
end


export getSSORCGDoubleTVRegularizationPreconditioner,wDoubleTVNodal,applySSORCGDoubleTV

function getSSORCGDoubleTVRegularizationPreconditioner(C,M::AbstractMesh,Iact, omega::Float64=1.0,tol::Float64=1e-2,maxCGIter::Int64=100)
	return HessianPreconditioner(SSORCGParamFourthOrder(C,M,Iact,omega,tol,maxCGIter),applySSORCGDoubleTV,setupDoNothing);
end

function wDoubleTVNodal(m::Vector, mref::Vector, M::AbstractMesh; Iact=1.0, C=[])	
	dm = m.-mref;
	d2R = wTVRegNodal(m,mref,M,Iact = Iact,C = C)[3];
	clear(M);
	d2R = d2R'*d2R;
	dR  = d2R*dm;
	Rc  = 0.5*dot(dm,dR);
   return Rc,dR,d2R
end

function applySSORCGDoubleTV(Hs::Function, d2R::SparseMatrixCSC,v::Vector,param)
	d2R = wTVRegNodal(v,v,param.M,Iact = param.Iact, C = param.C)[3]
	Dinv = param.omega./diag(d2R);
	aux   = zeros(size(d2R,2));
	SSOR(r) = (aux[:]=0.0; ssorPrecTrans!(d2R,aux,r,Dinv); return aux);
	x = KrylovMethods.cg(d2R,v,tol=param.tol,maxIter=param.maxCGIter,M=SSOR,out=-1)[1]
	x = KrylovMethods.cg(d2R',copy(x),tol=param.tol,maxIter=param.maxCGIter,M=SSOR,out=-1)[1]
	return x;
end





export getSensTMatVec

function getSensTMatVec(v::Vector,m::Vector,pFor::FWIparam)

    # extract pointers
    M  				= pFor.Mesh
    omega 			= pFor.omega
    gamma 			= pFor.gamma
    Q     			= pFor.Sources
    P     			= pFor.Receivers
    Ainv			= pFor.ForwardSolver
	
	batchSize 		= pFor.forwardSolveBatchSize;
	select  		= pFor.sourceSelection;
	
	nsrc = size(Q,2); nrec = size(P,2);
	
	if length(select) > 0
		nsrc = length(select);
	end
	
	if batchSize > nsrc
		batchSize = nsrc;
	end
	
	numBatches 	= ceil(Int64,nsrc/batchSize);
	n = prod(M.n+1);
	
	# ALL AT ONCE CODE
	
	H = spzeros(Complex128,n,n);
	if isa(Ainv,ShiftedLaplacianMultigridSolver)
		H = GetHelmholtzOperator(M,m,omega, gamma, true,true);
		Ainv = updateParam(Ainv,M,m,omega);
		H = H + GetHelmholtzShiftOP(m, omega,Ainv.shift); 
		H = H';
		# H is actually shifted laplacian now...
	elseif isa(Ainv,JuliaSolver)
		H = GetHelmholtzOperator(M,m,omega, gamma, true,true);
	end
	
	JTv = zeros(Float64,n)
	Vdatashape = reshape(v,nrec,nsrc);
	
	for k_batch = 1:numBatches
		batchIdxs = (k_batch-1)*batchSize + 1 : min(k_batch*batchSize,nsrc);
		V = P*Vdatashape[:,batchIdxs];
		V,Ainv = solveLinearSystem(H,V,Ainv,1);
		if useFilesForFields
			filename = getFieldsFileName(omega);
			file     = matopen(filename);
			U        = read(file,string("Ubatch_",k_batch));
			close(file);
		else
			U = pFor.Fields[:,batchIdxs];
		end 
		
		V = conj((1+1im*vec(gamma)).*U).*V; 
		
		# for jj=1:length(batchIdxs)
			# jj_b = batchIdxs[jj];
			# for ii = 1:n
				# @inbounds V[ii,jj] *= conj((1.0+1im*gamma[ii])*U[ii,jj]); # *U[ii,jj_b]
			# end
		# end
		JTv += omega^2*vec(real(sum(V,2)));
	end
	
	if isa(Ainv,ShiftedLaplacianMultigridSolver)
		# println("Clearing!!!");
		clear!(Ainv.MG);
	end
	
	# V = P*Vdatashape;
	# V,Ainv = solveLinearSystem(H,V,Ainv,1)   # Lam = ForwardSolver\(P*V);
		
	# # JTv    = conj((1+1im*vec(gamma)).*U).*V; 
	# # V      = 0;
	# # JTv	 = omega^2*vec(sum(real(JTv),2));

	# for jj=1:size(U,2)		
		# for ii = 1:size(U,1)
			# @inbounds V[ii,jj] *= conj((1.0+1im*gamma[ii])*U[ii,jj]);
		# end
	# end
	# JTv	  = omega^2*vec(real(sum(V,2)));
    return JTv
end


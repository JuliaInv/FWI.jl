export getData

function getData(m,pFor::FWIparam,doClear::Bool=false)

    # extract pointers
    M       	= pFor.Mesh
    omega   	= pFor.omega
	wavelet 	= pFor.WaveletCoef;
    gamma   	= pFor.gamma
    Q       	= pFor.Sources
    P       	= pFor.Receivers
	Ainv    	= pFor.ForwardSolver;
	batchSize 	= pFor.forwardSolveBatchSize;
	select  	= pFor.sourceSelection;
    
    nrec  		= size(P,2) 
    nsrc  		= size(Q,2)
	
	
	# allocate space for data and fields
	n = prod(M.n+1);
	# ALL AT ONCE DIRECT CODE
	H = GetHelmholtzOperator(M,m,omega, gamma, true,useSommerfeldBC);
	
	if isa(Ainv,ShiftedLaplacianMultigridSolver)
		Ainv = updateParam(Ainv,M,m,omega);
		H = H + GetHelmholtzShiftOP(m, omega,Ainv.shift); 
		H = H';
		# H is actually shifted laplacian now...
		Ainv.MG.relativeTol *= 1e-4;
	end
   
	
	if select==[]
		Qs = Q*wavelet;
	else
		Qs = Q[:,select]*wavelet;
	end
	
	nsrc 		= size(Qs,2);
	
	if batchSize > nsrc
		batchSize = nsrc;
	end
	
	
	Fields = [];
	
	if doClear==false
		if pFor.useFilesForFields 
			tfilename = getFieldsFileName(omega);
			tfile     = matopen(tfilename, "w");
		else
			Fields    = zeros(FieldsType,n   ,nsrc);
		end
	end
	
	numBatches 	= ceil(Int64,nsrc/batchSize);
	D 			= zeros(Complex128,nrec,nsrc);
	U 			= zeros(Complex128,n,batchSize);
	
	Ainv.doClear = 1;
	for k_batch = 1:numBatches
		batchIdxs = (k_batch-1)*batchSize + 1 : min(k_batch*batchSize,nsrc);
		if length(length(batchIdxs))==batchSize
			U[:] = full(Qs[:,batchIdxs]);
		else
			U = convert(Array{Complex128},full(Qs[:,batchIdxs]));
		end
		U,Ainv = solveLinearSystem(H,U,Ainv,0)
		
		
		# pic = round(Int64,rand()*1000);
		# u = reshape(real(U[:,1]),M.n[1]+1,M.n[2]+1);
		# figure(pic)
		# imshow(u'); colorbar();
		
		Ainv.doClear = 0;
		D[:,batchIdxs]      = (P'*U);
		
		if doClear==false
			if pFor.useFilesForFields 
				write(tfile,string("Ubatch_",k_batch),convert(Array{Complex64},U));
			else
				Fields[:,batchIdxs] = U;
			end
		end
	end
	
	# pic = round(Int64,rand()*1000);
	# figure(pic)
	# imshow(real(D')); colorbar();
	
	if isa(Ainv,ShiftedLaplacianMultigridSolver)
		Ainv.MG.relativeTol *= 1e+4;
	end
	
	pFor.ForwardSolver = Ainv;
	
	if doClear==false
		if pFor.useFilesForFields 
			close(tfile);
		else
			pFor.Fields = Fields;
		end
	end
	
	if doClear
		clear!(pFor);
	elseif isa(Ainv,ShiftedLaplacianMultigridSolver)
		clear!(Ainv.MG); 
	end
    return D,pFor
end


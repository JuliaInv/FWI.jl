export freqCont;

function freqCont(mc, pInv::InverseParam, pMis::Array{RemoteChannel},contDiv::Array{Int64}, windowSize::Int64,
			resultsFilename::String,dumpFun::Function,mode::String="Joint",startFrom::Int64 = 1,cycle::Int64=0,method::String="projGN")
Dc = 0;
for freqIdx = startFrom:(length(contDiv)-1)
	if mode=="1stInit"
		reqIdx1 = freqIdx;
		if freqIdx > 1
			reqIdx1 = max(2,freqIdx-windowSize+1);
		end
		reqIdx2 = freqIdx;
	elseif mode=="Joint"
		reqIdx1 = freqIdx;
		if freqIdx > 1
			reqIdx1 = max(1,freqIdx-windowSize+1);
		end
		reqIdx2 = freqIdx;
	end
	currentProblems = contDiv[reqIdx1]:contDiv[reqIdx2+1]-1;
	println("\n======= New Continuation Stage: selecting batches: ",reqIdx1," to ",reqIdx2,"=======\n");
	pMisTemp = pMis[currentProblems];
	pInv.mref = mc[:];	
	## SETTING UP STOCHASTIC SELECTION
	
	# if freqIdx == 1 && length(currentProblems)>1
		# for ii = currentProblems[2:end]
			# setSourceSelectionRatio(pForTemp[ii], 0.2);
			
		# end
	# else
		# for ii = currentProblems
			# setSourceSelectionRatio(pForTemp[ii], 0.2);
		# end
	# end
	Temp = splitext(resultsFilename);
	if cycle==0
		Temp = string(Temp[1],"_FC",freqIdx,"_GN",Temp[2]);
	else
		Temp = string(Temp[1],"_Cyc",cycle,"_FC",freqIdx,"_GN",Temp[2]);
	end
	function dumpGN(mc,Dc,iter,pInv,PF)
		dumpFun(mc,Dc,iter,pInv,PF,Temp);
	end
	Dc = 0;
	if method == "projGN"
		mc,Dc,flag = projGNCG(mc,pInv,pMisTemp,dumpResults = dumpGN);
	elseif method == "barrierGN"
		mc,Dc,flag = barrierGNCG(mc,pInv,pMisTemp,rho=1.0,dumpResults = dumpGN);
	end
	clear!(pMisTemp);
end
return mc,Dc;
end




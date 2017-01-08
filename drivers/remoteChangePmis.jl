function updateWd(pMis::Array{RemoteChannel},Dc::Array{RemoteChannel})
@sync begin
	@async begin
		for k=1:length(pMis)
			pMis[k] = remotecall_fetch(updateWd,pMis[k].where,pMis[k],Dc[k]);
		end
	end
end
return pMis;
end

function updateWd(pMisRF::RemoteChannel,Dc::RemoteChannel)
pMis  = take!(pMisRF)
Dc = fetch(Dc);
pMis.Wd = 1./(real(Dc - pMis.dobs) + 1e-3*mean(abs(pMis.dobs[:]))) + 1im./(imag(Dc - pMis.dobs) + 1e-3*mean(abs(pMis.dobs[:])));
put!(pMisRF,pMis)
return pMisRF;
end


function multWd(pMis::Array{RemoteChannel},beta::Float64)
@sync begin
	@async begin
		for k=1:length(pMis)
			pMis[k] = remotecall_fetch(multWd,pMis[k].where,pMis[k],beta);
		end
	end
end
return pMis;
end

function multWd(pMisRF::RemoteChannel,beta::Float64)
pMis  = take!(pMisRF)
pMis.Wd *= beta;
put!(pMisRF,pMis)
return pMisRF;
end
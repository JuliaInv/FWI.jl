function updateWd(pMis::Array{RemoteRef{Channel{Any}}},Dc::Array{RemoteRef{Channel{Any}}})
@sync begin
	@async begin
		for k=1:length(pMis)
			pMis[k] = remotecall_fetch(pMis[k].where,updateWd,pMis[k],Dc[k]);
		end
	end
end
return pMis;
end

function updateWd(pMisRF::RemoteRef{Channel{Any}},Dc::RemoteRef{Channel{Any}})
pMis  = take!(pMisRF)
Dc = fetch(Dc);
pMis.Wd = 1./(real(Dc - pMis.dobs) + 1e-3*mean(abs(pMis.dobs[:]))) + 1im./(imag(Dc - pMis.dobs) + 1e-3*mean(abs(pMis.dobs[:])));
put!(pMisRF,pMis)
return pMisRF;
end


function multWd(pMis::Array{RemoteRef{Channel{Any}}},beta::Float64)
@sync begin
	@async begin
		for k=1:length(pMis)
			pMis[k] = remotecall_fetch(pMis[k].where,multWd,pMis[k],beta);
		end
	end
end
return pMis;
end

function multWd(pMisRF::RemoteRef{Channel{Any}},beta::Float64)
pMis  = take!(pMisRF)
pMis.Wd *= beta;
put!(pMisRF,pMis)
return pMisRF;
end
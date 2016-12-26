module FWI
	
using jInv.Mesh
using jInv.Utils
using jInv.LinearSolvers
using jInv.InverseSolve
using EikonalInv
using KrylovMethods
using Multigrid
using ForwardHelmholtz
using MAT

import jInv.ForwardShare.getData
import jInv.ForwardShare.getSensTMatVec
import jInv.ForwardShare.getSensMatVec
import jInv.LinearSolvers.copySolver
	
import jInv.ForwardShare.ForwardProbType

FieldsType = Complex32

useSommerfeldBC = false;

fieldsFilenamePrefix = "tempFWIfields"

function getFieldsFileName(omega::Float64)
	omRound = string(round((omega/(2*pi))*100.0)/100.0);
	tfilename = string(fieldsFilenamePrefix,"_f",omRound,"_worker",myid(),".mat");
end


export clear!
export setSourceSelection,setSourceSelectionRatio,setSourceSelectionNum

export FWIparam, getFWIparam
type FWIparam <: ForwardProbType
    omega					:: Float64     # frequencies
	WaveletCoef				:: Complex128
    gamma					:: Vector{Float64}     # attenuation
    Sources					:: Union{Vector{Float64},SparseMatrixCSC,Array{Float64,2}}   # Sources
    Receivers				:: Union{Vector{Float64},SparseMatrixCSC,Array{Float64,2}}
	Fields					:: Array{FieldsType}
    Mesh      				:: RegularMesh
	ForwardSolver			:: AbstractSolver
	forwardSolveBatchSize	:: Int64
	sourceSelection			:: Array{Int64,1}
	useFilesForFields		:: Bool
end

function getFWIparam(omega::Float64, WaveletCoef::Complex128, gamma::Vector{Float64},
							Sources::Union{Vector{Float64},SparseMatrixCSC,Array{Float64,2}},
							Receivers::Union{Vector{Float64},SparseMatrixCSC,Array{Float64,2}},
							Mesh::RegularMesh, ForwardSolver:: AbstractSolver, workerList::Array{Int64},forwardSolveBatchSize::Int64=size(Sources,2),useFilesForFields::Bool = false)		
	return getFWIparam([omega], [WaveletCoef],gamma,Sources,Receivers, Mesh,ForwardSolver, workerList,forwardSolveBatchSize,useFilesForFields);
end

function getFWIparam(omega::Array{Float64}, WaveletCoef::Array{Complex128},gamma::Vector{Float64},
							Sources::Union{Vector{Float64},SparseMatrixCSC,Array{Float64,2}},
							Receivers::Union{Vector{Float64},SparseMatrixCSC,Array{Float64,2}},
							Mesh::RegularMesh, ForwardSolver::AbstractSolver, workerList::Array{Int64},forwardSolveBatchSize::Int64=size(Sources,2),useFilesForFields::Bool = false)
	
	continuationDivision = zeros(Int64,length(omega)+1);
	continuationDivision[1] = 1;
	
	if workerList==[]
		ActualWorkers = workers();
	else
		ActualWorkers = intersect(workerList,workers());
		if length(ActualWorkers)<length(workerList)
			warn("FWI: workerList included indices of non-existing workers.")
		end
	end
	numWorkers = length(ActualWorkers);
	pFor   = Array(RemoteChannel,numWorkers*length(omega));
	SourcesSubInd = Array(Array{Int64,1},numWorkers*length(omega));
	for k=1:length(omega)
		getFWIparamInternal(omega[k],WaveletCoef[k], gamma,Sources,Receivers,zeros(FieldsType,0), Mesh, 
									ForwardSolver, forwardSolveBatchSize ,ActualWorkers,pFor,(k-1)*numWorkers+1,SourcesSubInd,useFilesForFields);
		continuationDivision[k+1] = k*numWorkers+1;
	end
	return pFor,continuationDivision,SourcesSubInd # Array of Remote Refs
end




function getFWIparamInternal(omega::Float64, WaveletCoef::Complex128,gamma::Vector{Float64},
							Sources::Union{Vector{Float64},SparseMatrixCSC,Array{Float64,2}},
							Receivers::Union{Vector{Float64},SparseMatrixCSC,Array{Float64,2}},
							fields::Array{FieldsType}, Mesh::RegularMesh,
							ForwardSolver:: AbstractSolver, forwardSolveBatchSize::Int64,
							Workers::Array{Int64}, pFor::Array{RemoteChannel},startPF::Int64,SourcesSubInd::Array{Array{Int64,1},1},useFilesForFields::Bool = false)
	i = startPF; nextidx() = (idx=i; i+=1; idx)
	nsrc  = size(Sources,2);
	numWorkers = length(Workers);
	# send out jobs
	@sync begin
		for p=Workers
			@async begin
				while true
					idx = nextidx()
					if idx > startPF + numWorkers - 1
						break
					end
					I_k = getSourcesIndicesOfKthWorker(numWorkers,idx - startPF + 1,nsrc)
					SourcesSubInd[idx] = I_k;
					pFor[idx] = initRemoteChannel(getFWIparamInternal,p, omega,WaveletCoef,  gamma, Sources[:,I_k], Receivers, fields, Mesh, 
																			copySolver(ForwardSolver),forwardSolveBatchSize,useFilesForFields);
					wait(pFor[idx]);
				end
			end
		end
	end
	return pFor # Array of Remote Refs
end

function getFWIparamInternal(omega::Float64,WaveletCoef::Complex128,
							gamma::Vector{Float64},
	                        Sources::Union{Vector{Float64},SparseMatrixCSC,Array{Float64,2}},
							Receivers::Union{Vector{Float64},SparseMatrixCSC,Array{Float64,2}},
							Fields::Array{FieldsType},
							Mesh::RegularMesh, ForwardSolver:: AbstractSolver, forwardSolveBatchSize::Int64,useFilesForFields::Bool = false)
	return FWIparam(omega,WaveletCoef,gamma,Sources,Receivers,Fields,Mesh,ForwardSolver,forwardSolveBatchSize,Array(Int64,0),useFilesForFields)
end

# function setSourceSelection(pForRF::RemoteRef{Channel{Any}}, selection::Array{Int64,1})
	# s = copy(selection);
	# if minimum(s) < 1
		# error("FWI: sources selection out of range");
	# end
	# pFor  = take!(pForRF);
	# s = [];
	# if isa(pFor,FWIparam)
		# Sources = pFor.Sources
		# if maximum(s)>size(Sources,2)
			# s = s[s.<=size(Sources,2)];
			# warn("FWI: reducing selection: s = s[s.<=size(Sources,2)]");
		# end	
		# pFor.sourceSelection = s;
	# end
	# put!(pForRF,pFor);
	# return s;
# end


# function setSourceSelectionRatio(pForRF::RemoteRef{Channel{Any}}, selectionRatio::Float64)
	# if selectionRatio <= 0.0 || selectionRatio > 1.0
		# error("selectionRatio has to be between 0 and 1");
	# end
	# if selectionRatio == 1.0 
		# return;
	# end
	# pFor  = take!(pForRF);
	# s = [];
	# if isa(pFor,FWIparam)
		# Q = pFor.Sources
		# s = randperm(size(Q,2))[1:ceil(Int64,size(Q,2)*selectionRatio)];
		# pFor.sourceSelection = s;
	# end
	# put!(pForRF,pFor);
	# return s;
# end

# function setSourceSelectionNum(pForRF::RemoteRef{Channel{Any}}, selectionNum::Int64)
	# if selectionNum <= 0 
		# error("selectionNum has to be bigger than 1");
	# end
	# pFor  = take!(pForRF);
	# s = [];
	# if isa(pFor,FWIparam)
		# Q = pFor.Sources
		# s = randperm(size(Q,2))[1:selectionNum];
		# pFor.sourceSelection = s;
	# end
	# put!(pForRF,pFor);
	# return s;
# end


import jInv.Utils.clear!
function clear!(pFor::FWIparam)
	clear!(pFor.ForwardSolver);
	pFor.Fields = zeros(FieldsType,0);
	clear!(pFor.Mesh);
	return pFor;
end

include("getData.jl")
include("getSensMatVec.jl")
include("getSensTMatVec.jl")
include("FourthOrderHesPrec.jl")
include("freqCont.jl")
#include("timeDomainFWI.jl")
end
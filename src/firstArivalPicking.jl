export getFirstArival
function getFirstArival(dataTraces::Array{Array{Float32,2},1},mask,x_rcv_loc::Array{Float64},x_src_loc::Array{Float64},vmax::AbstractFloat,vmin::AbstractFloat,window::Int64,dt::Float64,offset)

dobsEikPicked = zeros(Float64,size(mask));

for k = 1:length(dataTraces)
	Dk = dataTraces[k];
	sigma = maximum(abs(Dk[:]));
	# print("Mean Signal: ")
	# println(sigma)
	
	dobsEik_k = vec(dobsEikPicked[:,k]);
	Dk = copy(Dk);

	Dk[:,3:end-2] = (Dk[:,1:end-4] + 1.5*Dk[:,2:end-3] + 2*Dk[:,3:end-2] + 1.5*Dk[:,4:end-1] + Dk[:,5:end] )./7;
	Dk[:,3:end-2] = (Dk[:,1:end-4] + 1.5*Dk[:,2:end-3] + 2*Dk[:,3:end-2] + 1.5*Dk[:,4:end-1] + Dk[:,5:end] )./7;
	dobsEikPicked[:,k] = pickTime(Dk,mask[:,k],x_rcv_loc,x_src_loc[k],vmax,vmin,window,dt,offset); # Try to compute travel time
end

return dobsEikPicked;
end





function pickTime(dataTrace::Array,mask::Array,x_rcv_loc::Array{Float64},x_src_loc::Float64,vmax::AbstractFloat,vmin::AbstractFloat,window::Int64,dt::Float64,zeroTimeOffset::Int64)
dpredEik = zeros(Int64,size(dataTrace,1));
println("window = ",window);
println("zeroTimeOffset = ",zeroTimeOffset);
println("vmax = ",vmax);
println("vmin = ",vmin);
window = window;
halfWindow = div(window,2);

for ii = 1:length(dpredEik)
	
	if mask[ii]==0 || x_rcv_loc[ii] < x_src_loc
		continue;
	end
	windowIdx = ii;
	lookAtPrev = ii > 1 && mask[ii-1]!=0;
	lookAtNext = ii < length(dpredEik) &&  mask[ii+1]!=0;
	if lookAtPrev && lookAtNext
		windowIdx = (ii-1):(ii+1);
	elseif lookAtPrev
		windowIdx = (ii-1):ii;
	elseif lookAtNext
		windowIdx = ii:(ii+1);
	end
	
	jjmin = 3*window;
	jjmax = size(dataTrace,2)-3*window;
	determineZeroTimeOffset = 0;
	if lookAtPrev && x_rcv_loc[ii-1] >= x_src_loc
		# jjmin = dpredEik[ii-1] + convert(Int64,round((x_rcv_loc[ii] - x_rcv_loc[ii-1])/(vmax*dt)));
		# jjmax = dpredEik[ii-1] + convert(Int64,round((x_rcv_loc[ii] - x_rcv_loc[ii-1])/(vmin*dt)));
	else
		# jjmin = convert(Int64,round((x_rcv_loc[ii] - x_src_loc)/(vmax*dt)))+zeroTimeOffset;
		determineZeroTimeOffset=1;
		# we know that this point is very close to the source and we know the velocity in the sea.
		# jjmax = convert(Int64,round((x_rcv_loc[ii] - x_src_loc)/(vmin*dt)))+zeroTimeOffset;
		# jjmin = jjmax;
	end


	vals = zeros(max(jjmax - jjmin +1,0));
	valsDeriv = zeros(max(jjmax - jjmin +1,0));
	# println(size(s_tag))

	for jj = jjmin:jjmax
		nearWindowSignal = dataTrace[windowIdx,(jj-window+1):jj];
		farWindowSignal = dataTrace[windowIdx,(jj-2*window):(jj-window)];
		vals[jj - jjmin + 1 ] = vecnorm(nearWindowSignal)/vecnorm(farWindowSignal);
		
		forwardWindowSignal = sum(vec(dataTrace[windowIdx,(jj+1):jj+halfWindow+1]));
		backwardWindowSignal = sum(vec(dataTrace[windowIdx,(jj-halfWindow):jj]));
		valsDeriv[jj - jjmin + 1] = forwardWindowSignal - backwardWindowSignal;
		
		
		# println("placing vals: ")
		if jj > jjmin
			if vals[jj - jjmin + 1 ] > 2.0 && valsDeriv[jj - jjmin + 1]*valsDeriv[jj - jjmin] < 0.0
				dpredEik[ii] = jj;
				break;
			end
		end
		
	end
	
	if determineZeroTimeOffset==1
		zeroTimeOffset = dpredEik[ii] - convert(Int64,round((x_rcv_loc[ii] - x_src_loc)/(vmin*dt)));
	end
	if dpredEik[ii] == 0
		warn("Could not find travel time - maybe recorded time is too short..., Aborting");
		break;
	end
	
end

for ii = length(dpredEik):-1:1
	if mask[ii]==0 || x_rcv_loc[ii] >= x_src_loc
		continue;
	end
	windowIdx = ii;
	lookAtPrev = ii > 1 && mask[ii-1]!=0;
	lookAtNext = ii < length(dpredEik) &&  mask[ii+1]!=0;
	if lookAtPrev && lookAtNext
		windowIdx = (ii-1):(ii+1);
	elseif lookAtPrev
		windowIdx = (ii-1):ii;
	elseif lookAtNext
		windowIdx = ii:(ii+1);
	end
	
	jjmin = 3*window;
	jjmax = size(dataTrace,2)-3*window;
	determineZeroTimeOffset = 0;
	if lookAtNext && x_rcv_loc[ii+1] < x_src_loc
	else
		determineZeroTimeOffset=1;
	end
	


	vals = zeros(max(jjmax - jjmin +1,0));
	valsDeriv = zeros(max(jjmax - jjmin +1,0));
	for jj = jjmin:jjmax
		nearWindowSignal = dataTrace[windowIdx,(jj-window+1):jj];
		farWindowSignal = dataTrace[windowIdx,(jj-2*window):(jj-window)];
		vals[jj - jjmin + 1 ] = vecnorm(nearWindowSignal)/vecnorm(farWindowSignal);
		
		forwardWindowSignal = sum(vec(dataTrace[windowIdx,(jj+1):jj+halfWindow+1]));
		backwardWindowSignal = sum(vec(dataTrace[windowIdx,(jj-halfWindow):jj]));
		valsDeriv[jj - jjmin + 1] = forwardWindowSignal - backwardWindowSignal;
		
		
		# println("placing vals: ")
		if jj > jjmin
			if vals[jj - jjmin + 1 ] > 2.0 && valsDeriv[jj - jjmin + 1]*valsDeriv[jj - jjmin] < 0.0
				dpredEik[ii] = jj;
				break;
			end
		end
	end
	
	if determineZeroTimeOffset==1
		zeroTimeOffset = dpredEik[ii] - convert(Int64,round(-(x_rcv_loc[ii] - x_src_loc)/(vmin*dt)));
	end
	if dpredEik[ii] == 0
		warn("Could not find travel time - maybe recorded time is too short..., Aborting");
		break;
	end

end



dpredEik = convert(Array{Float64},dpredEik);
for ii=1:length(dpredEik)
	dpredEik[ii] = max(dpredEik[ii] - zeroTimeOffset,0)*dt;
end

return dpredEik
end

function pickTimeSimple(dataTrace::Array,mask::Array,x_rcv_loc::Array{Float64},x_src_loc::Float64,vmax::AbstractFloat,vmin::AbstractFloat,window::Int64,dt::Float64,offset::Int64)
dpredEik = zeros(Int64,size(dataTrace,1));
window = 2*window;
for ii = 1:length(dpredEik)
	if mask[ii]==0 
		continue;
	end
	windowIdx = ii;
	lookAtPrev = ii > 1 && mask[ii-1]!=0;
	lookAtNext = ii < length(dpredEik) &&  mask[ii+1]!=0;
	if lookAtPrev && lookAtNext
		windowIdx = (ii-1):(ii+1);
	elseif lookAtPrev
		windowIdx = (ii-1):ii;
	elseif lookAtNext
		windowIdx = ii:(ii+1);
	end
	jjmin = 2*window+1;
	jjmax = size(dataTrace,2) - 4*window

	vals = zeros(jjmax - jjmin + 1);
	for jj = jjmin : jjmax
		windowSignal = vec(dataTrace[windowIdx,jj-2*window:(jj-window)]);
		wideWindowSignal = vec(dataTrace[windowIdx,(jj-window+1):jj]);
		stdShort = norm(windowSignal);
		stdLong = norm(wideWindowSignal);
		vals[jj - jjmin + 1 ] = abs(stdLong)./abs(stdShort);
		if vals[jj - jjmin + 1 ] > 10.0
			dpredEik[ii] = jj;
			found = 1;
			break;
		end
	end
	dpredEik[ii] = jjmin + indmax(vals) + 1;
end
dpredEik = convert(Array{Float64},dpredEik);
for ii=1:length(dpredEik)
	dpredEik[ii] = max(dpredEik[ii] - offset,0)*dt;
end

return dpredEik
end





############ OLD IDEA FOR BACKUP

# for ii = length(dpredEik):-1:1
	# if mask[ii]==0 || x_rcv_loc[ii] >= x_src_loc
		# continue;
	# end
	# windowIdx = ii;
	# lookAtPrev = ii > 1 && mask[ii-1]!=0;
	# lookAtNext = ii < length(dpredEik) &&  mask[ii+1]!=0;
	# if lookAtPrev && lookAtNext
		# windowIdx = (ii-1):(ii+1);
	# elseif lookAtPrev
		# windowIdx = (ii-1):ii;
	# elseif lookAtNext
		# windowIdx = ii:(ii+1);
	# end
	
	# jjmin = 1;
	# jjmax = size(dataTrace,2);
	# if lookAtNext && x_rcv_loc[ii+1] < x_src_loc
		# jjmin = dpredEik[ii+1] + convert(Int64,round((x_rcv_loc[ii+1] - x_rcv_loc[ii])/(vmax*dt)));
		# jjmax = dpredEik[ii+1] + convert(Int64,round((x_rcv_loc[ii+1] - x_rcv_loc[ii])/(vmin*dt)));
	# else
		### jjmin = convert(Int64,round(-(x_rcv_loc[ii] - x_src_loc)/(vmax*dt)))+zeroTimeOffset-10;
		# jjmax = convert(Int64,round(-(x_rcv_loc[ii] - x_src_loc)/(vmin*dt)))+zeroTimeOffset;
		# jjmin = jjmax;
	# end
	# found = 0;
	# jjmax = min(jjmax,size(dataTrace,2));
	# found = 0;
	# vals = zeros(max(jjmax - jjmin +1,0));
	# for jj = jjmin:jjmax
		# farWindowSignal = vec(dataTrace[windowIdx,(jj-2*window):(jj-window)]);
		# nearWindowSignal = vec(dataTrace[windowIdx,(jj-window+1):jj]);
		# stdNear = norm(nearWindowSignal);
		# stdFar = norm(farWindowSignal);
		# vals[jj - jjmin + 1 ] = (abs(stdNear)./abs(stdFar));
		# if vals[jj - jjmin + 1 ] > 4
			# dpredEik[ii] = jj;
			# found = 1;
			# break;
		# end
	# end
	### println(vals);
	# if found == 0 && jjmax >= jjmin
		# dpredEik[ii] = jjmin + indmax(vals) + 1;
	# end
	
	# if dpredEik[ii] == 0
		# warn("Could not find travel time - maybe recorded time is too short..., Aborting");
		# break;
	# end

# end



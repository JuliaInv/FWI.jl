function getAnalyticalConstGrad2D(n::Array{Int64,1},h::Array{Float64,1})
src = [1,div(n[2],2)];

source1 = (src[1]-1)*h[1];
source2 = (src[2]-1)*h[2];

(X1,X2) = ndgrid((0:(n[1]-1))*h[1],(0:(n[2]-1))*h[2]);

a  = -0.4;
s0 = 2.0;

kappaSquared = s0^2 + 2*a*(X1-source1);

SBarSquared = s0^2 + a*(X1-source1);
rSquared = ((X1-source1).^2 + (X2-source2).^2);
sigmaSquarred = 2*rSquared./(SBarSquared + sqrt(SBarSquared.^2 - (a^2).*rSquared));
sigma = sqrt(sigmaSquarred);
T_exact =  SBarSquared.*sigma - (a^2).*(sigma.^3)./6; 

return kappaSquared,src,T_exact;
end

function getAnalyticalConstGrad3D(n::Array{Int64,1},h::Array{Float64,1})
src = [div(n[1],2),div(n[2],2),1];

source1 = (src[1]-1)*h[1];
source2 = (src[2]-1)*h[2];
source3 = (src[3]-1)*h[3];

(X1,X2,X3) = ndgrid((0:(n[1]-1))*h[1],(0:(n[2]-1))*h[2],(0:(n[3]-1))*h[3]);

# a  = -1.75;
a  = -1.65;
s0 = 2.0;

kappaSquared = s0^2 + 2*a*(X3-source3);
SBarSquared = s0^2 + a*(X3-source3);
rSquared = (X1-source1).^2 + (X2-source2).^2 + (X3-source3).^2;
sigmaSquarred = 2*rSquared./(SBarSquared + sqrt(SBarSquared.^2 - (a^2).*rSquared));
sigma = sqrt(sigmaSquarred);
T_exact =  SBarSquared.*sigma - (a^2).*(sigma.^3)./6; 

return kappaSquared,src,T_exact;
end



function runAnalyticalConstGradExperiment3D(h::Array{Float64,1},WU::Float64)
I = [1.6,1.6,0.8];

n = zeros(Int64,3);
n[1] = round(Int64,(I[1]/h[1])+1);
n[2] = round(Int64,(I[2]/h[2])+1);
n[3] = round(Int64,(I[3]/h[3])+1);

# src = [div(n[1],2),div(n[2],2),div(n[3],2)];
src = [div(n[1],2),div(n[2],2),1];

source1 = (src[1]-1)*h[1];
source2 = (src[2]-1)*h[2];
source3 = (src[3]-1)*h[3];

(X1,X2,X3) = ndgrid((0:(n[1]-1))*h[1],(0:(n[2]-1))*h[2],(0:(n[3]-1))*h[3]);

a  = -1.75;
s0 = 2.0;

kappaSquared = s0^2 + 2*a*(X3-source3);
SBarSquared = s0^2 + a*(X3-source3);
rSquared = (X1-source1).^2 + (X2-source2).^2 + (X3-source3).^2;
sigmaSquarred = 2*rSquared./(SBarSquared + sqrt(SBarSquared.^2 - (a^2).*rSquared));
sigma = sqrt(sigmaSquarred);
T_exact =  SBarSquared.*sigma - (a^2).*(sigma.^3)./6; 

sigmaSquarred = sigma = SBarSquared = X1 = X2 = X3 = rSquared = 0.0;

return T_exact,src,kappaSquared;
end




function getAnalyticalConstGradInv2D(n::Array{Int64,1},h::Array{Float64,1})

src = [1,div(n[2],2)];

source1 = (src[1]-1)*h[1];
source2 = (src[2]-1)*h[2];

(X1,X2) = ndgrid((0:(n[1]-1))*h[1],(0:(n[2]-1))*h[2]);

a  = 1.0;
s0 = 2.0;

Kappa = 1./(1./s0 + a*(X1-source1));
T_exact = (1/a)*acosh(1+(0.5*s0*a*a).*Kappa.*((X1-source1).^2 + (X2-source2).^2));
G2_exact = (1/a)*(1./sqrt((1+(0.5*s0*a*a).*Kappa.*((X1-source1).^2 + (X2-source2).^2)).^2 - 1)).*(0.5*s0*a*a).*Kappa.*2.*(X2-source2);
G1_exact = (1/a)*(1./sqrt((1+(0.5*s0*a*a).*Kappa.*((X1-source1).^2 + (X2-source2).^2)).^2 - 1)).*((0.5*s0*a*a).*Kappa.*2.*(X1-source1) 
			+ (0.5*s0*a*a).*((X1-source1).^2 + (X2-source2).^2).*(-(Kappa.^2)).*a);
kappaSquared = Kappa.^2;

return kappaSquared,src,T_exact;

end


function getAnalyticalConstGradInv3D(n::Array{Int64,1},h::Array{Float64,1})

src = [div(n[1],2),div(n[1],2),1];

source1 = (src[1]-1)*h[1];
source2 = (src[2]-1)*h[2];
source3 = (src[3]-1)*h[3];

(X1,X2,X3) = ndgrid((0:(n[1]-1))*h[1],(0:(n[2]-1))*h[2],(0:(n[3]-1))*h[3]);

a  = 1.0;
s0 = 2.0;
Kappa = 1./(1./s0 + a*(X3-source3));
Rsquared = (X1-source1).^2 + (X2-source2).^2 + (X3-source3).^2;
T_exact = (1/a)*acosh(1+(0.5*s0*a*a).*Kappa.*(Rsquared));
G1_exact = (1/a)*(1./sqrt((1+(0.5*s0*a*a).*Kappa.*Rsquared).^2 - 1)).*(0.5*s0*a*a).*Kappa.*2.*(X1-source1);
G2_exact = (1/a)*(1./sqrt((1+(0.5*s0*a*a).*Kappa.*Rsquared).^2 - 1)).*(0.5*s0*a*a).*Kappa.*2.*(X2-source2);
G3_exact = (1/a)*(1./sqrt((1+(0.5*s0*a*a).*Kappa.*Rsquared).^2 - 1)).*((0.5*s0*a*a).*Kappa.*2.*(X3-source3) 
			+ (0.5*s0*a*a).*(Rsquared).*(-(Kappa.^2)).*a);
kappaSquared = Kappa.^2;

return kappaSquared,src,T_exact;

end


function getSmoothGaussianMedium(n::Array{Int64,1},h::Array{Float64,1})
src = div(n,4);
kappaSquared = [];
T_exact = [];
if length(n)==2
	(T1_exact,G11_exact,G12_exact) = getSmoothFactoredModel(n,h);
	(T0,G01,G02,L0) = getAnalytic2DeikonalSolutionAll(n,h,src);
	G1_exact = T0.*G11_exact + G01.*T1_exact;
	G2_exact = T0.*G12_exact + G02.*T1_exact;
	kappaSquared = G1_exact.*G1_exact + G2_exact.*G2_exact;
	T_exact = T0.*T1_exact;
else

	(T1_exact,G11_exact,G12_exact,G13_exact) = getSmoothFactoredModel3D(n,h);

	(T0,G01,G02,G03) = getAnalytic3DeikonalSolutionAll(n,h,src);

	G1_exact = T0.*G11_exact + G01.*T1_exact;
	G2_exact = T0.*G12_exact + G02.*T1_exact;
	G3_exact = T0.*G13_exact + G03.*T1_exact;
	kappaSquared = G1_exact.*G1_exact + G2_exact.*G2_exact + G3_exact.*G3_exact;
	T_exact = T0.*T1_exact;
	
end

return kappaSquared,src,T_exact;
end

function getSmoothFactoredModel(n::Array{Int64,1},h::Array{Float64,1})
xsrc = zeros(2);
src_kappa = zeros(Int64,2);
src_kappa[1] = div(n[1],3);
src_kappa[2] = div(n[2],4);

xsrc[1] = (src_kappa[1]-1)*h[1];
xsrc[2] = (src_kappa[2]-1)*h[2];

X1,X2 = ndgrid((0:(n[1]-1))*h[1],(0:(n[2]-1))*h[2]);

sigma = 0.1;
T1_exact = (exp( - (sigma*((X1 - xsrc[1]).^2) + 4*sigma*((X2-xsrc[2]).^2))) + 1)/2;
G11_exact = -2*sigma*(X1 - xsrc[1]).*exp( - (sigma*((X1 - xsrc[1]).^2) + 4*sigma*((X2-xsrc[2]).^2)))/2;
G12_exact = -8*sigma*(X2 - xsrc[2]).*exp( - (sigma*((X1 - xsrc[1]).^2) + 4*sigma*((X2-xsrc[2]).^2)))/2;

return T1_exact,G11_exact,G12_exact;
end

function getSmoothFactoredModel3D(n::Array{Int64,1},h::Array{Float64,1})
xsrc = zeros(3);
src_kappa = zeros(Int64,3);
src_kappa[1] = div(n[1],3);
src_kappa[2] = div(n[2],4);
src_kappa[3] = div(n[3],2);

xsrc[1] = (src_kappa[1]-1)*h[1];
xsrc[2] = (src_kappa[2]-1)*h[2];
xsrc[3] = (src_kappa[3]-1)*h[3];

X1,X2,X3 = ndgrid((0:(n[1]-1))*h[1],(0:(n[2]-1))*h[2],(0:(n[3]-1))*h[3]);

Sigma = [0.1,0.4,0.2];
EXPRSkewedSquared = exp(-Sigma[1]*((X1 - xsrc[1]).^2) - Sigma[2]*((X2-xsrc[2]).^2) - Sigma[3]*((X3-xsrc[3]).^2));

T1_exact = EXPRSkewedSquared/2 + 1/2;
G11_exact = -Sigma[1]*(X1 - xsrc[1]).*EXPRSkewedSquared;
G12_exact = -Sigma[2]*(X2 - xsrc[2]).*EXPRSkewedSquared;
G13_exact = -Sigma[3]*(X3 - xsrc[3]).*EXPRSkewedSquared;

return T1_exact,G11_exact,G12_exact,G13_exact;
end

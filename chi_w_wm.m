%%%% Calculating the real and Matsubara frequency spin susceptibility
%%%% chi%%%%%

a=3.86;                                  %lattice constant (A)
hbar=6.58211928*10^(-16);
kB=8.6173303*10^-5;                     %Boltzmann constant (eV/K)
T=70;                                   %adjust for corresponding data (K)
beta=1/(kB*T);                          %1/eV
                          
                        
%We first start in real frequencies, where we use the RPA to calcualte chi. 
%To this end, we first need to load bare suseptibility (chi0_T) in real frequencies
%first, calculated in the script realw_chi0.

U=0.660;     %adjust for corresponding data set

chi0=chi0_70;       %adjust for each temperature
dims=size(chi0);
N=dims(1)*dims(2);
om=linspace(0.0,0.4,dims(3));
dw=om(2)-om(1);
qx_=linspace(-1/(2*a),1/(2*a),dims(2));   qy_=linspace(-1/(2*a),1/(2*a),dims(2));   %(1/m)
[qx,qy] = meshgrid(qx_,qy_); 
Uq=-U*(cos(2*pi*qx*a)+cos(2*pi*qy*a))/2;       %t-J - like coupling function

chi=chi0./(1-bsxfun(@times,Uq,chi0));

chi=imag(chi);     %keep only the imaginary part, which is the one measured in INS
                   %and entering the Lehmann representation

chi=reshape(chi,dims(1)*dims(2),dims(3));
schi=zeros(N,dims(3));

for i=1:N
	schi(i,:)=-fliplr(chi(i,:));
end

chi=horzcat(schi,chi);chi(:,dims(3))=[];
chi=reshape(chi,dims(1),dims(2),2*dims(3)-1);

%%%%%%%%%%%%%%% -------------------%%%%%%%%%%%%%%%%%
% chi is the interacting spin susceptibility, sampled over a 128 x 128
%momentum grid and covering 320 energies from -0.4 to 0.4. It is odd in
%energy. In the paper we plot this quantity at pi,pi for various
%temperatures
%%%%%%%%%%%%%%% -------------------%%%%%%%%%%%%%%%%%

om1=horzcat(-fliplr(om),om);om1(dims(3))=[];

%Now calculate using Lehmann representation

wmmax=3.49;                      %eV ; This is enough to capture the decay in Chi 
mmax=wmmax*beta/pi; mmax=round(mmax);
if rem(mmax,2)~=0,
    mmax=mmax+1;
end

m=-mmax:2:mmax;
wm=m*pi/beta;              %Bosonic Matsubara frequencies (eV)
Nm=length(m);

for i=1:Nm
    for j=1:length(om1)
        integ(:,:,j)=-(1/pi)*chi(:,:,j)/(1i*wm(i)-om1(j));
    end
    chiM(:,:,i)=dw*trapz(integ,3);                  %units of sec
end

%These two steps are NOT compulsory, but may be needed when the number of
%bosonic Matsubara frequencies computed here does not match the number of
%Fermionic frequencies computed for the Green's function in script
%"matsubara_G_BSE"
chiM(:,:,end)=[];chiM(:,:,1)=[];
wm(end)=[];wm(1)=[];


%For the Lehmann representation, the wm=0 points will have divergences because the integrand is 0/0 when
%om=0. Therefore, we calculate chiM(k,wm=0) separately and then re-introduce
%them in chiM. This is easily done since Imchi is linear near om=0, so around
%this point the integrand becomes (a*om)/om=a, where a is the slope.

dims=size(chiM);

for i=1:length(om1)
    integ(:,:,i)=(1/pi)*chi(:,:,i)/om1(i);
end

slope(:,:)=chi(:,:,floor((length(om1)/2)+1)+1)/dw;   %slopes of Im chi near w=0 at all momenta
integ(:,:,floor((length(om1)/2)+1))=(1/pi).*slope;

chiM(:,:,floor(dims(3)/2)+1)=dw*trapz(integ,3);
chiM=real(chiM);

%%%%%%   chiM is the Matsubara spin susceptibility, sampled over a 128 x 128 
%momentum grid and with Nm bosonic Matsubara energies (a number which depends 
%on temperature), spanning energies from -3.5 to 3.5 eV   %%%%%%%






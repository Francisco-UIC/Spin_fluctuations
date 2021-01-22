%This script essentially reproduces the results of Chatterjee, et. al. on
%the Dynamic spin susceptibility in Bi2212. The calculations are thus performed
%in REAL frequency.

tic

a=3.86;                               %lattice constant (A)
kB=8.6173303*10^-5;                   %Boltzmann constant (eV/K)
T=70;                                 %adjust for corresponding data (K)
beta=1/(kB*T);                        %1/eV
fbin=323;

%ARPES Data: the data have been symmetrized (to kill the Fermi function)
%and a k-independent background spectrum has been subtracted. Moreover, 
%they have been to sammple 1 quadrant of the BZ evenly.

%bscesT variable should be loaded in workspace first !!!!! 
%(T is the termperature)

A=eval(strcat('bsces',num2str(T)));             %spectral function as function of kx, ky, and f (3 dimension)
dims=size(A);

e=(1:dims(3))-fbin;
dw=0.002;               %data were sampled at 0.002 eV energy intervals
e=e*dw;                 %binding energy scale

%%%%% ----First, we reduce momentum and energy sampling density of the
%%%%%               spectral function -------           %%%%%%%

%The present calculation can be quite time-consuming, so I de-interpolate
%the data by producing values of A on a coarser k-space grid.

[X,Y]=meshgrid(1:dims(1));
ki=65;                  %reduced # of momentum points in Brillouin Zone quadrant
[Xq,Yq]=meshgrid(linspace(1,dims(1),ki));
Ai=zeros(ki,ki,dims(3));

for i=1:dims(3)
    Ai(:,:,i)=interp2(X,Y,A(:,:,i),Xq,Yq);
end

A=Ai;           %new spectral function with reduced momentum sampling density
clear Ai 
dims=size(A);


%We also don't need such dense energy sampling, so I also de-interpolate
%the energies to a fraction of the measured sampling density.

frac=2;
edc=reshape(A,dims(1)*dims(2),dims(3));
ei=round(dims(3)/frac);
wi=linspace(E(1),E(end),ei);
edci=zeros(dims(1)*dims(2),ei);

for i=1:dims(1)*dims(2)
    edci(i,:)=interp1(E,edc(i,:),wi);
end

Ai=reshape(edci,dims(1),dims(2),ei);

A=Ai;           %new spectral function with reduced energy sampling density
dw=dw*frac;

%Perform symmetry operations on quadrant to cover whole BZ

sA=flip(A,2);             
sA(:,end,:)=[];           
sA2=horzcat(sA,A);        
ssA2=flip(sA2,1);         
ssA2(end,:,:)=[];         
A=vertcat(ssA2,sA2); 
dims=size(A);
N=dims(1)*dims(2);

%The energy-momentum sum of the spectral function must equal N:

tsw=dw*sum(sum(trapz(A,3)))/N;
A=A/tsw;

clear sA sA2 ssA2

%%%% ------ Now, we can proceed to calculate particle-hole bubble ------%%%%

f=1./(exp(wi*beta)+1);            %Fermi function
d=0.0055;       %eV  convergence factor in chi0

corr=zeros(dims(1),dims(2),dims(3));
integ=zeros(dims(1),dims(2),dims(3));
sum1=zeros(dims(1),dims(2),dims(3)); 

om=linspace(0,0.4,160);     %160 energies ("OMega") over which chi0 will be sampled (eV)
om=0;
chi0=zeros(dims(1),dims(2),length(om));

%The momentum sum in chi0 is a convolution, so we handle it with Fourier
%techniques. Thereore, we start by calculating the 2-D momentum F-transform
%of the spectral functions:
Arw=fftshift(fft2(ifftshift(A)));

%Now, for the sums:

for k=1:length(om)          %energies of chi0
    for i=1:dims(3)         %first energy dummy variable
        for j=1:dims(3)     %second energy dummy variable
            corr(:,:,j)=fftshift(ifft2(ifftshift(Arw(:,:,i).*Arw(:,:,j))));    %convolution theorem for momentum sums
            integ(:,:,j)=corr(:,:,j)*(f(i)-f(j))/(om(k)+wi(i)-wi(j)+1i*d);     %integrand for energy sums
        end
        sum1(:,:,i)=dw*sum(integ,3);    %first energy sum
    end
    chi0(:,:,k)=-(dw/N)*sum(sum1,3);    %second energy sum. 1/N falls out of the fft2 definition
end

chi0_70=chi0;

%%%% chi0 is the bare susceptibility (particle-hole bubble), sampled over
%%%% the BZ and covering 160 energies in the range [0,400] meV %%%%%

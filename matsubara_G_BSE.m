%%%% Matusbara Green;s function and solution of the Bethe-Salpeter Equation %%%%


%First, we must construct the full spectral function using its symmetry.
%This is done in another script (ces.m); the data set must be
%ready before using this script.

%some useful constants:

a=3.86;                              %lattice constant (A)
kB=8.6173303*10^-5;                  %Boltzmann constant (eV/K)
T=70;                                %adjust for corresponding data (K)
beta=1/(kB*T);                       %1/eV
fbin=323;                            %"bin" corresponding to Fermi energy for thin film data (UD67 and OP80)
%fbin=261;                           %"bin" corresponding to Fermi energy for single crystal data (OP91)

%ARPES Data: the data have been symmetrized (to kill the Fermi function)
%and a k-independent background spectrum has been subtracted. Moreover, 
%they have been to sammple one quadrant of the BZ evenly.

%bscesT variable should be loaded in workspace first !!!!! 
%(T is the termperature) 

A=eval(strcat('bsces',num2str(T)));            %spectral function as function of kx, ky, and f (3 dimension)
dims=size(A);

e=(1:dims(3))-fbin;
dE=0.002;           %data were sampled at 0.002 eV energy intervals
e=e*dE;  
E=single(e);

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
clear Ai Xq Yq  
dims=size(A);

N=dims(1)*dims(2);

%We also don't need such dense energy sampling, so I also de-interpolate
%the energies to a fraction of the measured sampling density.

edc=reshape(A,dims(1)*dims(2),dims(3));
frac=2;
ei=round(dims(3)/frac);
wi=linspace(E(1),E(end),ei);
edci=zeros(dims(1)*dims(2),ei);

for i=1:dims(1)*dims(2)
    edci(i,:)=interp1(E,edc(i,:),wi);
end

Ai=reshape(edci,dims(1),dims(2),ei);
dims=size(Ai);

A=Ai;       %new spectral function with reduced energy sampling density
dE=dE*frac;

%Normalize the TOTAL spectral weight to unity:

tsw=dE*sum(sum(trapz(A,3)))/N;
A=A/tsw;

%Set fermionic frequency range

wpmax=15;           %eV (G drops below 5% of its max)

pmax=wpmax*beta/pi; pmax=round(pmax);
if rem(pmax,2)==0,
    pmax=pmax+1;
end

p=-pmax:2:pmax;
wp=p*pi/beta;                  %fermionic Matsubara frequencies (eV

Np=length(wp);

integ=zeros(dims(1),dims(2),dims(3));
G=zeros(dims(1),dims(2),length(wp));

%With the spectral function, we calculate the Green's function via the 
%Lehmann representation. For the calculations ahead, we focused on the
%normal state, so no anomalous Green's functions were necessary
 
for i=1:Np
    for j=1:length(wi)
        integ(:,:,j)=A(:,:,j)/(1i*wp(i)-wi(j));
    end
    G(:,:,i)=dE*trapz(integ,3);                  %units of sec
end

clear evec sevec integ 

G=1i*imag(G);             
%this follows from symmetrization of the spectral function  A(k,w)=A(k,-w)
%near the Fermi surface

%Now, get G for the entire BZ.

sG=flip(G,2);             
sG(:,end,:)=[];          
sG2=horzcat(sG,G);        
ssG2=flip(sG2,1);         
ssG2(end,:,:)=[];         
Gkw=vertcat(ssG2,sG2);    

clear sG sF sG2 sF2 ssG2 ssF2
dims=size(Gkw);
N=dims(1)*dims(2);


% Now, solve the Bethe-Salpeter Equation using the power method

P=-Gkw.*Gkw;      %"Pairing Kernel" (again, follows from symmetrization)

%This pairing kernel decays very fast (wp^-2), so we introduce an
%energy cutoff at 3.5 eV. We have checked that cutoffs at larger energies
%make no difference on the final results for P and the eigenvalues.

wp=wp-wp(Np/2+1);
[val highcut]=min(abs(wp-3.5));
[val lowcut]=min(abs(wp+3.5));
wp(highcut:end)=[];        wp(1:lowcut)=[];
Np=length(wp);

P(:,:,highcut:end)=[];  P(:,:,1:lowcut)=[];

n=20;        %number of iterations for power method

%Construct starting guess Phi_0:

kx_=linspace(-1/(2*a),1/(2*a),dims(1));     %(1/m)
[kx,ky] = meshgrid(kx_,kx_); 
clear kx_ 
phi0=cos(2*pi*kx*a)-cos(2*pi*ky*a); 

%starting guess for eigenfunction (d-wave symmetry, with no energy dependence)
%see Monthoux, Phd Thesis

lambda=zeros(1,n);

%Again, we use Fourier transforms to compute convolutions
%Start with the full 3-D F. transform of the Matsubara chi.
chirt=fftshift(fftn(ifftshift(chiM)));
chirt=real(chirt);                
%chiM are the Matusubara spin susceptibilities obtained in script
%"chi_w_wm"at each temperature.  LOAD THEM FIRST!!!!
%By definition, these are positive, real, and even in 
%both momenta and Matsubara energy. Hence, all F-transforms are real, but 
%numerics may introduce small imaginary parts.

phi{1}=phi0;
lambda(1)=1;            %starting guess for eigenvalue
U=0.66;                %spin-fermion coupling; obtained from script "chi_w_wm", at each T.
for i=1:n
    prod=bsxfun(@times,P,phi{i});            %form a product which is to be convolved with chiM
    Prt=fftshift(fftn(ifftshift(prod))); 
    phi{i+1}=-1.5*bsxfun(@times,U^2,(1/(N*beta))*fftshift(ifftn(ifftshift(chirt.*Prt))));    %BSE, solved with 3-D Fourier convolution theoerm
    ali1=max(max(max(abs(phi{i+1}))));    ali=max(max(max(abs(phi{i}))));
    lambda(i+1)=(ali+ali1)/2;            %average j-th and (j+1)-th eigenvalues, for stability (not strictly necessary)
    phi{i+1}=phi{i+1}/lambda(i+1); 
    phi{i+1}=real(phi{i+1});             %phi comes out real, but picks up small artificial imaginary part 
    phi{i}=[];        %don't store eigenfunctions at each iteration, only keep last one 
end

%%%% phi is the pairing eigenfunction, sampled over a 128 x 128 momentum
%%%% grid (covering the entire BZ) and with Np Matsubara points (a number
%%%% which depends on temperature) spanning energies from -3.5 to 3.5 eV %%%%%

%%%% Also, lambda for the n+1-th iteration is the leading eigenvalue of the
%%%% Bethe Salpeter Equation %%%%%

plot(1:n+1,lambda(:),'color',rand(1,3));  %plot convergence behavior and read-off last eigenvalue
ylabel('Leading eigenavalue'); xlabel('Iteration number');

%%
%Check desired eigenfunction symmetries

iter=n+1;
ebin=ceil(Np/2);     %eV (this is the "wn=0" bin depends on Temperature) 
figure;

%Extract phi along (pi,0 - 0,pi) diagonal (at zeroth Matsubara frequency) and plot it
phik=real(phi{iter}(:,:,ebin));
phik(ceil(dims(1)/2):end,:)=[];  phik(:,ceil(dims(2)/2):end)=[];
phik=rot90(phik);
phid=diag(phik);  
plot(1:length(phid),phid,'color',rand(1,3));
ylabel('\Phi along (\pi,0)-(0,\pi) direction'); 

%Extract phi(pi,0,w_n) and plot it
phiw=real(phi{iter}(1,65,:)); phiw=reshape(phiw,1,Np);
figure;
plot(wp,phiw,'color',rand(1,3));
ylabel('\Phi(\pi,0,\omega_n)'); xlabel('\omega_n');


   
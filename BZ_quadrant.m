%%%% Preparing the ARPES data for analysis  %%%%

% This script is representative of the data handling that was required
% before the calculations of the spin response and the BSE could be
% performed. This script accomplishes various tasks:

% 1) Puts the ARPES data in the form of a 3-D array, two
% dimensions for momenta, the third for binding energy. These data are
% assumed to be loaded in the workspace. Moreover, they have already been stored
% in a cell array consisting of 19 cells, where each cell represents the
% dispersion along the ky direction for a given kx direction.

% 2) Interpolates the data to effectively sample the data over an evenly
% spaced (with equal density along kx and ky directions) BZ.

% 3) Subtracts background and superlattice signals, normalizes, and 
% symmetrizes all spectra.


%%%% This script is specific to the UD80 sample at 70K; slightly different
%%%% paramters were sometimes used for different samples.

%%%%%%%%%%%%%%%% ------------------------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ef = 17.644;          % (eV) Fermi energy determined from Au spectrum. Au was in 
                      % electrical contact with Bi2212 sample, so they have                      % the same chemical potential
fbin = 323;           % "bin" corresponding to Fermi energy

data70k = data;                 % "data" are the data already loaded in the workspace
dims = size(data70k);

% Transform into regular 3-d array, for ease of manipulations

data70k = cell2mat(data70k); 
data70k = reshape(data70k,esteps,ksteps,dims(2));
data70k = permute(data70k,[2,3,1]);

% Eliminate 1st wave as it contains Au information (see data log).

data70k(:,1,:) = [];
dims = size(data70k);

% The following two manipulations re-orient the data in momentum space, for
% convenience.


for j = 1:dims(3)
    data70k(:,:,j) = fliplr(data70k(:,:,j));
    data70k(:,:,j) = flipud(data70k(:,:,j));
end
%%
% Now, I will subtract the background. To avoid subtracting a noisy spectrum
% from an already noisy data set, I will pick a subset of points near the BZ 
% corner, where the signal is mostly background, then I'll average the spectra  
% and subtract this from all data, since the background is almost
% momentum independent.

subdat = data70k(1:30,1:2,:);      
dimsub = size(subdat);
bg = mean(subdat,1);
bg = reshape(bg,dimsub(2),dimsub(3));
bg = mean(bg,1);
bg = bg-min(bg);

% Before subtracting background from each spectrum, normalize it to every
% spectrum:

edc = reshape(data70k,dims(1)*dims(2),dims(3));
r = bg(1)./edc(:,1);
edc = bsxfun(@times,edc,r);

% Now that we have the background, we subtract it from all spectra:

for i = 1:dims(1)*dims(2)
    edc(i,:) = edc(i,:)-bg;
end

dE = e(2)-e(1);

data70k = reshape(edc,dims(1),dims(2),dims(3));

% After subtraction, there is still some superlattice contribution and some
% remaining noisy signal throughout the BZ, which we eliminate as follows. 

for i = 1:dims(3)
    ces = data70k(:,:,i);
    ces(ces<0.2*max(max(ces))) = 0;    
    data70k(:,:,i) = ces; 
end
clear ces
%%
% Next, I interpolate the data to have equal sampling density on kx and ky:

kx = round(dims(2)*dims(1)/(k(end)-k(1)));   % k variable was created in "import"; it contains the photoemission angles
fsi = zeros(dims(1),kx,dims(3));
x = 1:dims(2);
xi = linspace(1,dims(2),kx);

for j = 1:dims(3) 
    for i = 1:dims(1) 
        fsi(i,:,j) = interp1(x,data70k(i,:,j),xi);
    end
end

% Since we want to isolate one quadrant of the BZ, and our data extends
% symmetrically over the X and Y quadrants of the BZ, we eliminate the -ky
% half of the data.

fsi(118:end,:,:) = [];        
dims = size(fsi);
fsi(isnan(fsi)==1) = [];

dims = size(fsi);

edc = reshape(fsi,dims(1)*dims(2),dims(3));
edc = fsi;

%%% Now, we have the data sampled through the Y quadrant of the zone.
%%% Since the spectral function is also symmetric with respect to the BZ
%%% diagonals, we remove the data above the kx=ky diagonal. To this end, we
%%% first find the SC node on the Fermi surface, construct the
%%% diagonal that contains it, and eliminate the data above this diagonal.

nodx = 54; nody = 122;                   %nodal position, ascertained by individually 
                                     %looking at EDCs on Fermi surface   
clear edc2 fsi

xp = abs(nody-nodx);
yp = abs(dims(1)-(dims(2)-xp));   %amount of padding in x and y directions

dims = size(edc);

for j = 1:dims(3) 
    
    %Now, cut off data beyond the node, i.e., above (nody-nodx)-th
    %diagonal.
    
    a0 = tril(edc(:,:,j),abs(nody-nodx));
    a1 = tril(edc(:,:,j),abs(nody-nodx)-1);
    
    % I'll use padding to extend the remaining arrays a and a1 to a dimension such
    % that the (nody-nodx)-th diagonal becomes the main diagonal 
   
    pad = padarray(flipud(a0),[xp,yp],'post');
    ces(:,:,j) = fliplr(pad);
end

clear  a a1 pad pad1 quad xp yp
dims = size(ces);

%%%%% "ces" are the Constant Energy Surfaces at KINETIC energies in the range
%%%%% [17,17.75] eV. Later on we'll subtract the Fermi energy to put this
%%%%% scale in terms of BINDING energies. These CESs are non-zero in an
%%%%% irreducible octant of the BZ; later on we'll reflect about the
%%%%% symmetry axes to cover the (almost) entire BZ.    

%%% -------------------------------------------------------------------------- %%%

%%% Symmetrization %%%

% This section contains the commands for proper symmetrization of the
% spectral function (see Chatterjee et al, PHYSICAL REVIEW B 75, 172504). To
% this end, a complicated geometrical operation is necessary, involving 
% momentum vectors perpendicular to the Fermi surface. To simplify the
% operation, I used a TB fit to the dispersion, and used the resulting FS as
% a reference.

% Compare with Norman's dispersion:

a = 3.86;

kx_ = linspace(0,1/(2*a),dims(1));   ky_ = linspace(0,1/(2*a),dims(1));

[kx,ky] = meshgrid(kx_,ky_); 

% For the bare band ek, we use Norman's parametrization.

mu = 19;   % sometimes a small shift is needed to make the TB fit match the FS

% tb2
c1 = 0.1960;  c2 = -0.6798;  c3 = 0.2368;  c4 = -0.0794;  c5 = -0.0343;  c6 = 0.0011;

ek = mu+1000*(c1+(c2/2)*(cos(2*pi*kx*a)+cos(2*pi*ky*a))+c3*(cos(2*pi*kx*a).*cos(2*pi*ky*a))...
    +(c4/2)*(cos(2*pi*2*kx*a)+cos(2*pi*2*ky*a))+(c5/2)*(cos(2*pi*2*kx*a).*cos(2*pi*ky*a)+...
    cos(2*pi*kx*a).*cos(2*pi*2*ky*a))+c6*(cos(2*pi*2*kx*a).*cos(2*pi*2*ky*a))); %meV
ek1 = ek;
ek1(ek1<-2)  =0;
ek1(ek1>2) = 0;     
ek1 = abs(triu(ek1));      % for working only with octant           
[row,col] = find(ek1);     % read out coordinates of tight binding fit FS
fsp = length(row);
ek2 = zeros(dims(1),dims(2));

for i = 1:fsp
    ek2(row(i),col(i)) = 1; 
end

mesh(ces(:,:,fbin)); view(2); hold on;  mesh(2*ek2);view(2)  % check that fit actually passes through FS


%%% Symmetrization implies that the spectral functionat any point k can be 
% written as A(k,w)=I(k,w)+I(2kf-k,-w), where kf is the Fermi momentum at
% which k is perpendicular to the Fermi surface. "I" is the ARPES data we
% have already reduced to the "ces" variable.

% So, first, we will need I(-w), which can easily get by flipping our
% EDCs from "left to right" (or w <--> -w):

dims = size(ces);

edc = reshape(ces,dims(1)*dims(2),dims(3));
edc = padarray(edc,[0 2*fbin-dims(3)-1],'post');
edc2 = fliplr(edc);
edc = reshape(edc,dims(1),dims(2),2*fbin-1);
edc2 = reshape(edc2,dims(1),dims(2),2*fbin-1);

dims = size(edc);

% Now we use the equation above, by noting that kf is really the point in the FS
% closest to the point of interest k.

d = 1000*ones(dims(1),dims(2));
bsces70 = zeros(dims(1),dims(2),dims(3));

for h = 1:dims(1)      % runs over kx points
    for j = 1:h        % runs along ky up to diagonal
        k1 = [j,h];    % create k vector for specific j,h
        for i = 1:fsp
            fsc = [row(i),col(i)];            % FS coordinates
            d(row(i),col(i)) = norm(k1-fsc);  % distance from k to all FS points
        end
        minimum = min(min(d));
        [x,y] = find(d==minimum); kf = [x(end),y(end)];
        k3 = 2*kf-k1; b = k3<=0; s = sum(b);
        if s >= 1
        bsces70(j,h,:) = edc(j,h,:);
        else 
        bsces70(j,h,:) = edc(j,h,:) + edc2(k3(1),k3(2),:);   % actual symmetrization step
        end
    end
end


% Throw away EDCs too far beyond FS, since they cannot be properly
% constructed by symmetrizing
ecut = 80;                    %energy cutoff (meV)
b = ek>ecut;        

for i = 1:dims(3)
    cut = bsces70(:,:,i);  cut = reshape(cut,dims(1),dims(2));
    cut(b) = 0;
    bsces70(:,:,i) = cut;
end

%%% Now, the data cover some momenta outside the Fermi surface; moreover,
%%% they extend from -640 to 640 meV in binding energy.

% Reflect about diagonal

for i = 1:dims(3)
    ref = bsces70(:,:,i)';
    bsces70(:,:,i) = bsces70(:,:,i)-diag(diag(bsces70(:,:,i)));
    bsces70(:,:,i) = bsces70(:,:,i)+ref;
end

% The symmetrized data now span the whole Y-quadrant   %%%%%

% To finalize, we normalize the EDCs to unit area, as is required from
% the basic properties of the spectral function

sedc = reshape(bsces70,dims(1)*dims(2),dims(3));
e1 = e-e(fbin);
evec = e1(1:fbin);
sevec = fliplr(evec); evec(end)=[];
E = horzcat(evec,-sevec);
dE = E(2)-E(1);
sedc = bsxfun(@rdivide,sedc,dE*trapz(sedc,2));    
sedc(isnan(sedc)) = 0;

bsces70 = reshape(sedc,dims(1),dims(2),2*fbin-1);
dims = size(bsces70);

% The spectral fucntion thus created has some outlying intensities at
% momentum points where it should be zero. (Note that the normalization above 
% can amplify signals which were originally small)

bsces70(1:56,1:56,:) = 0;
bsces70(72,72,319:332) = min(min(bsces70(:,:,319:332))); 

%%%%%%%%%%%%%% --------------------%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUMMARY:

%  "bsces70" is basically the spectral function of the OP80 thin film
%  at 70K. It covers one quadrant of the BZ and binding in the range [-640,80] meV

%%%%%%%%%%%%%%%---------------------%%%%%%%%%%%%%%%%%%%%%%


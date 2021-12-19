% Import and pre-normalize data from specified folder.

% User should have a folder with 16-19 files (depending on data set in
% question). Each of these files corresponds to a dispersion relation with
% kx fixed and varying ky. These files are pulled directly from the
% detector software, and contain useful information, such as maximum and
% minimum electron kinetic energies, maximum and minimum photoemission
% angles, step time, number of energy and momentum points, etc. This script
% extracts all these parameters and stores them in variables to be used in
% another script. Of course, this script also extracts the actual
% ARPES data and stores them in the form of a cell array with 16-19 cells.
% This output array is called "data".

% Retrieve files using 'dir' function

% Here, user should provide path to folder containing data files. The user 
% should make sure there are no files with same names under the same parent 
% directory, otherwise MATLAB might open the wrong file. Moreover, the
% folder must be in the current MATLAB directory.

folder = uigetdir;
files = dir(folder);
files(1:2)=[];            
names = { files.name };

% Get number of files

dims = size(names);
numfiles = dims(2);
clear dims;

% Open all files:

for i = 1:numfiles
    
    data{i} = dlmread(names{i},'',40,0);      %(40 for thin film data)
  
    % Get number of sweeps and divide by it:
    
    file = names{i};
    fid = fopen(file,'r');
    
    for j = 1:38                       % 38 for thin film data       
        numsweep{i} = fgetl(fid);
    end
    numsweep{i}(1:17) = [];
    numsweep{i} = str2num(numsweep{i});
    
    % Also, get step time, since it affects overall counts:
    
    fid = fopen(file,'r');
   
    for j = 1:30                %(30 for thin film data)
    stept{i} = fgetl(fid);
    end
    stept{i}(1:10) = [];
    stept{i} = str2num(stept{i});
    
    % First, subtract lowest value of tail:
    
    data{i} = bsxfun(@minus,data{i}(:,1:end),min(data{i},[],1));
     
    % Since each wave is taken at a different photon flux, they must each be
    % normalized independently to a value that does not change between
    % slices. Since the data are noisy, this edge value fluctuates.
    % Therefore, we average the spectra over a few energies close to this reference
    % energy.

    dims = size(data{i});
    subdat = data{i}(1:5,:);
    av = mean(subdat,1);
  
    data{i} = bsxfun(@rdivide,data{i}(:,1:end),av);
    data{i}(:,1) = [];           
   
end

% Extract relevant parameters, such as angular
% range, energy step size, etc. from raw data file

file = names{1};

% Find minimum and maximum energies and number of energy points

fid = fopen(file,'r');
for i = 1:9
    eline = fgetl(fid);
end
eline(1:18) = [];
e = str2num(eline); clear eline
emin = e(1);
emax = e(end);
edims = size(e);
esteps = edims(2);

% Find minimum and maximum momenta and number of momentum points

fid = fopen(file,'r');
for i = 1:12
    kline = fgetl(fid);
end
kline(1:18) = [];
k = str2num(kline);  clear kline
kmin = k(1);
kmax = k(end);
kdims = size(k);
ksteps = kdims(2);

% Find energy step size

fid = fopen(file,'r');

for i = 1:29                   % (29 for thin film data)
    deltaEline = fgetl(fid);
end
deltaEline(1:12) = [];
deltaE = str2num(deltaEline);

fclose(fid)

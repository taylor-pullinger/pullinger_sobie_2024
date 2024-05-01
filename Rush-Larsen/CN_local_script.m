clear
addpath('Rush-Larsen')
addpath('ToRORd19')
addpath('Helpers')
minerva = 0; % Allows adjusted settings for running on the Mount Sinai HPC, 
% "Minerva". To run on a local machine, set to 0.

% Choose population heterogeneity (our study mostly used 1 and 3) 
% 1) uniform cable of CMs from a population with sigma = 0.1, 
% 2) sigma = 0.05, 
% 3) sigma = 0.1, 
% 4) sigma = 0.2, 
% 5) sigma = 0.25, 
% 6) sigma = 0.3, 
% 7) sigma = 0.4
tissue = [3]; 

sample = 1; % select which cable(s) to run

% if applying a100 drug multiplier, what should it be? If not, set to 1
drugEffects = [1]; 

% which conductance is targeted by the drug? If no drug, set to 1
% order is stored in NamesCell.mat
drugTargets = [1];

% to create a uniform stretch in the center of a heterogeneous cable, 
% indicate it here, otherwise set to 50. 
% e.g. 47:53 will create a uniform stretch of cell #50 from position 
% 47 to 53. To retain the original cable length of 100 CMs, the original 
% CMs 1-3 & 98-100 are removed.
unistretch = 50;

% to order the CMs in a cable by a particular conductance set this variable
% to that conductance's number (positive for ascending, negative for
% descending). For example, to arrange CMs in descending order by GKr
% value, set ascend = -5.
ascend = 0;

% indicate which portion of the cable to apply the drug, as a 2-element
% vector, e.g. for upstream [1, 49], downstream [51, 100]. For the full
% cable set to [1, 100].
drug_loc = [1, 100];

% for parallelization efficiency, all conditions are run for the ic_s2c
% CM first (so slow collection of initial conditions can be done in parallel), 
% after which the further s2cs CMs may be run. Set as a fraction of cable length. 
ic_s2c =  50/100; 
s2cs = (46:55)/100;

% this script allows for looping through two variables simultaneously (e.g. 
% different drug strengths over a number of samples. Input those variables 
% (set above) in this statement.  
nloops = length(sample)*length(drugEffects); 

cores = nloops; % number of cores needed for parallel computation

% use this if-statement to limit the number of cores requested, for example
% to the number of cores on your machine
if cores > 8 
   cores = 8;
end

 % to run simulation(s) in series, comment out the next 3 lines
 delete(gcp('nocreate')) % delete any parallel pool active from a prior run
 pp = parpool(cores); % creates a new parallel pool with necessary cores
 
 % this option ensures 1 cable is calculated on 1 core
 opts = parforOptions(pp,'RangePartitionMethod','fixed','SubrangeSize',1);

parfor (i = 1:nloops, opts) % to run in series this line should read: for i = 1:nloops
   % index marker for iterated variable 1 (e.g. sample). Insert the name
   % of the other iterated variable (e.g. drugEffects) into this statement.
   si = ceil(i/length(drugEffects)); 

   % index marker for iterated variable 2 (e.g. drugEffects). Insert the name
   % of this iterated variable (e.g. drugEffects) into this statement.
   di = mod(i-1,length(drugEffects)) + 1;

   % this vector contains 4 elements:
   % 1) a toggle for whether this simulation will be calculating the VW, 1
   % for yes, 0 for no
   % 2) the S2-CM
   % 3) an initial guess at the time for VW opening in ms following S1
   % 4) an initial guess at the time for VW closing in ms following S1
   vw = [1, ic_s2c, 280, 291];%260, 310]; 

   drugs = struct(); % combine drug variables to pass to the function
   drugs.target = drugTargets(1);
   drugs.effect = drugEffects(di);
   drugs.unistretch = unistretch;
   drugs.drug_loc = drug_loc(1,:);

   % Rj multiplier & S2 amplitude as a 2-element vector
   rj_s2a = [5, 210];

   % pass all necessary data to the VW function.
   VW(vw, tissue(1), sample(si), minerva, drugs, ascend, rj_s2a)
end

% Check that the settings in the below loop are identical to the above to
% run additional S2-CMs. Comment out to take a single measurement.
%{ 
for j = 1:length(s2cs)
   parfor (i = 1:nloops, opts)
      si = ceil(i/length(drugEffects)); 
      di = mod(i-1,length(drugEffects)) + 1; 
      vw = [1, s2cs(j), 260, 310];
      drugs = struct();
      drugs.target = drugTargets(1);
      drugs.effect = drugEffects(di);
      drugs.unistretch=unistretch;
      drugs.drug_loc = drug_loc(1,:);
      rj_s2a = [5, 210];
      VW(vw, tissue(1), sample(si), minerva, drugs, ascend, rj_s2a) 
   end
end
%}

% if running in series, comment this line out
delete(pp)

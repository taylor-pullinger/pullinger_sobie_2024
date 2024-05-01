function [t, V_all, Cai_all, state_all_last_beat, ij_data] = f_tomek_cable_CNRL_CVTS(N, sample, tissue, protocol, vw)
% Simulation of propagating action potential using Tomek 2019 model
% Perform semi-implicit "Crank-Nicolson" solution
% Gating variables updated with Rush-Larsen method

if ~exist('vw','var')
   vw = [0 0 0 0];
end

if ~isfield(protocol, 'drugs')
   protocol.drugs = ones(N,17);
   protocol.drugs(protocol.drug_loc(1):protocol.drug_loc(2),protocol.drug_target) = protocol.drug_effect;
end

% ability to manipulate sodium gating-variables with drug targets 18, 19, &
% 20
tau_drugs = struct('m', 1, 'h', 1, 'j', 1);
if protocol.drug_target == 18
   tau_drugs.m = protocol.drug_effect;
elseif protocol.drug_target == 19
   tau_drugs.h = protocol.drug_effect;
elseif protocol.drug_target == 20
   tau_drugs.j = protocol.drug_effect;
end

if protocol.drug_target > 17
   protocol.drugs = ones(N,17) ;
end

if ~isfield(protocol, 'iso')
   protocol.iso = 0;
end

if ~isfield(protocol, 'downsample_step')
   protocol.downsample_step = 1;
end

if ~isfield(protocol, 'ascend')
   protocol.ascend = 0;
end
ascend = protocol.ascend;
if ascend > 0
   protocol.ascendtag = ['_ascendingG', num2str(ascend)];
elseif ascend < 0
   protocol.ascendtag = ['_descendingG', num2str(ascend*-1)];
else
   protocol.ascendtag = '';
end

if ~isfield(protocol, 'Rj_mod')
   protocol.Rj_mod = 1;
end
if protocol.Rj_mod ~= 5
   protocol.rjtag = ['_Rj_', num2str(protocol.Rj_mod)];
else
   protocol.rjtag = '';
end


if ~isfield(protocol, 'stim_cell')
   protocol.stim_cell = 1;
end

if ~isfield(protocol, 'stretch')
      stretchtag = '';
      protocol.stretch=[];
else
      stretchtag = protocol.stretchtag;
end

% Identify necessary file names
if tissue > 1
   tissue_string = num2str(tissue);
else
   % by default all uniform cables are created from a cell drawn from a 
   % population with sigma = 0.1. This can be changed by setting this 
   % variable, and its partner in VW.m (line 81) to
   % a different value
   uniform_variant = 3;
   tissue_string = num2str(tissue + uniform_variant * 0.1);
end
ic_filename = ['Rush-Larsen/ICs/ics_N', num2str(N),...
   '_t', tissue_string,...
   '_s', num2str(sample),...
   '_beat21', protocol.drugtag, protocol.ascendtag, protocol.rjtag, stretchtag,'.mat'];
pre_ic_filename = ['Rush-Larsen/ICs/preics_N', num2str(N),...
   '_t', tissue_string,...
   '_s', num2str(sample), protocol.drugtag, protocol.ascendtag, protocol.rjtag, stretchtag,'.mat'];

%% SETTINGS
n_stim = protocol.n_stim;

dt_small = protocol.dt_small; 
dt_big = protocol.dt_big; 

iso = protocol.iso;
drugs = protocol.drugs;

if vw(1)
   stim_cell = floor(vw(2) * N * 1.001);
   s2_time = protocol.s2_time;
   vw_ic_time = floor(s2_time);
   stim_delay = s2_time-vw_ic_time;
else
   stim_cell = protocol.stim_cell;
   stim_delay = protocol.stim_delay;
end

%%
init_Y = f_getInit_Y('endo');
% create a cardiomyocyte object using the ToRORd19 model (see
% ToRORd19_fast.mat)
CM = ToRORd19_fast(init_Y, 0, 0) ; 

%% Define constants that will not be varied randomly
%% Physical constants
F = CM.universals.F;  % Faraday constant, coulombs/mmol
R = CM.universals.R; % gas constant, J/K
T = CM.universals.T; % absolute temperature, K            

%% Cell Geometry
L = CM.geometry.L ;                                    % cm
rad = CM.geometry.rad ;                                % cm
Ageo = CM.geometry.Ageo ;                              % cm^2, geometric membrane area
Acap = CM.geometry.Acap ;                              % cm^2, capacitive membrane area

vcell = CM.geometry.vcell ;
vmyo = CM.geometry.vmyo ;
vnsr = CM.geometry.vnsr ;
vjsr = CM.geometry.vjsr ;
vss = CM.geometry.vss ;

%% Standard ionic concentrations for ventricular cells
ko = CM.universals.ko;
cao = CM.universals.cao;
nao = CM.universals.nao;
cli = 24;   % Intracellular Cl  [mM]
clo = 150;  % Extracellular Cl  [mM]

%% Cell constants
ICaL_fractionSS = 0.8 ;
INaCa_fractionSS = 0.35 ;

Cm = 1;
celltype = 0;

%% Set up pacing protocol
stim_dur = protocol.stim_dur ;
stim_amp = protocol.stim_amp ;
stim_freq = protocol.stim_freq ; % Hz

stim_times = zeros(2*n_stim + 1, 2);
stim_indices = zeros(2*n_stim + 1, 1);
stim_times(1,2) = stim_delay;
for i = 1:n_stim
   stim_times(2*i,2) = stim_times(2*i-1,2)+stim_dur;
   stim_times(2*i+1,2) = 1000/stim_freq*i+stim_delay;
   stim_times(2*i:2*i+1,1) = stim_times(2*i-1:2*i,2);
   stim_indices(2*i) = 1;
end

if vw(1)
   vw_run_dur = 80;
   load(ic_filename, 'state_all_s1', 't_s1')
   if max(squeeze(state_all_s1(:,1,end))) > -82
      vw_run_dur = 180;
   end
   [~, meanCV2080, ~] = calculateCV(t_s1', squeeze(state_all_s1(:,1,:))');
   if meanCV2080 < 12 % if propagation is slow, add time to capture sufficient data
      vw_run_dur = 120;

   end
   clear state_all_s1 t_s1 meanCV2080
   stim_times(end,2) = stim_times(end,1) + vw_run_dur ; 
end

tend = stim_times(end,2); % end of simulation, ms

slow_dur = protocol.slow_dur;
delay_chunk = 0 : dt_big : stim_delay * 1.000000001;
slow_chunk = stim_delay : dt_small : stim_delay* 1.000000001+slow_dur;
fast_chunk = stim_delay+slow_dur : dt_big : stim_delay* 1.000000001+1000/stim_freq;
chunk_lengths = [length(delay_chunk), length(slow_chunk)-1, length(fast_chunk)-1];
iterations = chunk_lengths(1) + n_stim * (chunk_lengths(2) + chunk_lengths(3));
t = delay_chunk;
for steps = 1:n_stim
   t = [t, (steps-1)*1000/stim_freq+[slow_chunk(2:end),fast_chunk(2:end)]];
end

if vw(1)
   t = t(t <= tend);
   iterations = length(t);
end

if dt_big > 1e-3 && length(t) < length(union(t,stim_times(:,1)))
   warning('pacing not captured precisely by time steps - results may be inaccurate')
end

if isfield(protocol, 'start_time')
   stim_times = [protocol.start_time, protocol.end_time];
   t = protocol.start_time:dt_small:protocol.end_time; 
   iterations = length(t);
   tend = stim_times(end,2);
end

%% Set up matrix for solving voltage equation implicitly
% Only has to be calculated once
% Cable parameters
% Variables defining cable properties

Rj = 1.5e3 * protocol.Rj_mod; % junctional resistance, kOhm (mV/kOhm = uA)

tda = -1/(Acap*Rj)*ones(N,1) ;
tdc = -1/(Acap*Rj)*ones(N,1) ;

tdb_slow = (2/(Acap*Rj) + Cm/dt_small)*ones(N,1) ;
tdb_slow(1) = 1/(Acap*Rj) + Cm/dt_small ;
tdb_slow(N) = 1/(Acap*Rj) + Cm/dt_small ;

A = spdiags([ tda tdb_slow tdc ], -1:1, N, N);
% makes a sparse NxN tridiagonal matrix:
% tdb(1) tdc(1) 0      0      0      0      0
% tda(2) tdb(2) tdc(2) 0      0      0      0
% 0      tda(3) tdb(3) tdc(3) 0      0      0
% 0      0      tda(4) tdb(4) tdc(4) 0      0
% 0      0      0      tda(5) tdb(5) tdc(5) 0
% 0      0      0      0      tda(6) tdb(6) tdc(6)

A_inv_slow = inv(A) ;
% inverts the matrix:
% A * A_inv = NxN identity matrix
% a matrix can be inverted when its determinant != 0 (matrix is linearly
% independent) & it is square

tdb_fast = (2/(Acap*Rj) + Cm/dt_big)*ones(N,1) ;
tdb_fast(1) = 1/(Acap*Rj) + Cm/dt_big ;
tdb_fast(N) = 1/(Acap*Rj) + Cm/dt_big ;

A = spdiags([ tda tdb_fast tdc ], -1:1, N, N);

A_inv_fast = inv(A) ;

%% Read in the population

load('Populations/minicombo.mat', 'sigmas')
sigma = sigmas(tissue);
if vw(1) % in vw case, just read in the relevent pre-specified population
   conductances = load(['Populations/sample',num2str(sample),'_vw_conds_c', num2str(N),'_s', num2str(sigma),'_t', tissue_string, stretchtag, '.mat'], 'cond_tensor');
   conductances = conductances.cond_tensor' .* drugs;
else
   conds_selection = load(['Populations/sample',num2str(sample),'_conds_selection.mat'], 'conds_selection');
   conds_selection = conds_selection.conds_selection;
   conductances = load('Populations/minicombo.mat','cond_tensor');
   if tissue > 1
      conductances = conductances.cond_tensor(conds_selection(1:N), :, tissue);
   else
      if sample > 1
         conductances = conductances.cond_tensor(conds_selection(1:N), :, uniform_variant);
      else
         conductances = conductances.cond_tensor(conds_selection(1:N), :, tissue);
      end
   end
   if tissue == 1
      conductances = repmat(conductances(1,:), N ,1); % make uniformly the first cell
   end

   if ~isempty(protocol.stretch)
      conductances = conductances';
      repCM=repmat(conductances(:,50),1,length(protocol.stretch));
      nToTrim=(length(protocol.stretch)-1)/2;
      upstreamHet = [conductances(:,nToTrim+1:49)];
      downstreamHet = [conductances(:,51:end-nToTrim)];
      conductances = [upstreamHet,repCM,downstreamHet];
      conductances = conductances';
   end

   conductances = conductances .* drugs;
end

if ascend
   conductances = sortrows(conductances, ascend);
end


%% Initial conditions
% Use the CN_getInits function to obtain initial conditions for the
% simulation
if isfile(ic_filename) % if we already have ICs
   if vw(1) % if doing s2 part of VW
      [v, nai, nass, ki, kss, cai, cass, cansr, cajsr, m, hp, h, j, jp, mL,...
         hL, hLp, a, iF, iS, ap, iFp, iSp, d, ff, fs, fcaf, fcas, jca, nca, nca_i,...
         ffp, fcafp, xs1, xs2, Jrel_np, CaMKt, ikr_c0, ikr_c1, ikr_c2, ikr_o, ikr_i,...
         Jrel_p] = CN_getInits(ic_filename, N, vw_ic_time, protocol);
      Istim = 0 * ones(N,1) ;
      last_beat_startdex = iterations + 1;
      state_all_last_beat = [];
      ij_data=[];
      t_final=t;

   elseif isfield(protocol, 'start_time')
         [v, nai, nass, ki, kss, cai, cass, cansr, cajsr, m, hp, h, j, jp, mL,...
         hL, hLp, a, iF, iS, ap, iFp, iSp, d, ff, fs, fcaf, fcas, jca, nca, nca_i,...
         ffp, fcafp, xs1, xs2, Jrel_np, CaMKt, ikr_c0, ikr_c1, ikr_c2, ikr_o, ikr_i,...
         Jrel_p] = CN_getInits(ic_filename, N, protocol.start_time, protocol);
      Istim = 0 * ones(N,1) ;
      last_beat_startdex = iterations + 1;
      state_all_last_beat = [];
      ij_data=[];
      t_final=protocol.end_time;


   else % if just doing a regular (non-VW) cable simulation
      [v, nai, nass, ki, kss, cai, cass, cansr, cajsr, m, hp, h, j, jp, mL,...
         hL, hLp, a, iF, iS, ap, iFp, iSp, d, ff, fs, fcaf, fcas, jca, nca, nca_i,...
         ffp, fcafp, xs1, xs2, Jrel_np, CaMKt, ikr_c0, ikr_c1, ikr_c2, ikr_o, ikr_i,...
         Jrel_p] = CN_getInits(ic_filename, N, [], protocol);
      Istim = 0 * ones(N,1) ;
      last_beat_startdex = find(t == stim_times(end-1,1));

      if n_stim > 5
         % for initial conditions, only save data at every ms timepoint to
         % save space
         downsampling = t(last_beat_startdex):protocol.downsample_step:t(end);
         [~,ia,~]=intersect(t,downsampling);
         if length(downsampling)~=length(ia)
            warning('downsampling error')
         end
         t_final=t(ia);
         state_all_last_beat = zeros(N,43,length(ia));
         ij_data = zeros(N,20,length(ia));
      else
         t_final = t;
         state_all_last_beat = [];
         ij_data=[];
      end
   end
else % if we don't have ICs yet
   save(pre_ic_filename, 'tissue_string', 'sample', 'tissue')
   protocol.tissue_string = tissue_string;
   protocol.sample = sample;
   protocol.tissue = tissue;
   [v, nai, nass, ki, kss, cai, cass, cansr, cajsr, m, hp, h, j, jp, mL,...
      hL, hLp, a, iF, iS, ap, iFp, iSp, d, ff, fs, fcaf, fcas, jca, nca, nca_i,...
      ffp, fcafp, xs1, xs2, Jrel_np, CaMKt, ikr_c0, ikr_c1, ikr_c2, ikr_o, ikr_i,...
      Jrel_p] = CN_getInits(pre_ic_filename, N, [], protocol);
   Istim = 0 * ones(N,1) ;
   last_beat_startdex = find(t == stim_times(end-1,1));
   downsampling = t(last_beat_startdex):t(end);
   [~,ia,~]=intersect(t,downsampling);
   if length(downsampling)~=length(ia)
      warning('downsampling error')
   end
   t_final=t(ia);
   state_all_last_beat = zeros(N,43,length(ia));
   ij_data = zeros(N,20,length(ia));
end

%% Make vectors for data to keep
if n_stim <= 5
   V_all = zeros(iterations,N) ;
   Cai_all = zeros(iterations,N) ;
else
   V_all = [] ;
   Cai_all = [] ;
end

%% set model constants outside loop
% CaMK constants
KmCaMK = 0.15 ;
aCaMK = 0.05 ;
bCaMK = 0.00068 ;
CaMKo = 0.05 ;
KmCaM = 0.0015 ;
% INab
PNab = 1.9239e-09 .* conductances(:,11) ;
% ICab
PCab = 5.9194e-08 .* conductances(:,12) ;
% IpCa
GpCa = 5e-04 .* conductances(:,13) ;
% chloride
ecl = (R .* T ./ F) .* log(cli ./ clo) ;            % [mV]

Fjunc = 1 ;% fraction in SS and in myoplasm - as per literature, I(Ca)Cl is in junctional subspace
Fsl = 1 - Fjunc ;

GClCa = conductances(:,14) .* 0.2843 ;   % [mS./uF]
GClB = conductances(:,15) .* 1.98e-3 ;        % [mS./uF]
KdClCa = 0.1 ;    % [mM]
% calcium buffer constants
cmdnmax = 0.05 ;
if celltype == 1
   cmdnmax = cmdnmax .* 1.3 ;
end
kmcmdn = 0.00238 ;
trpnmax = 0.07 ;

BSRmax = 0.047 ;
KmBSR = 0.00087 ;
BSLmax = 1.124 ;
KmBSL = 0.0087 ;
csqnmax = 10.0 ;
kmcsqn = 0.8 ;

%% Loop to run simulation
for its=1:iterations
   if iso
      cass = min(cass, 0.025 * ones(N,1));
   end

   t_now = t(its) ;
   if t_now < tend
      dt = t(its+1)-t(its);
   else
      dt = dt_big;
   end

   if mod(t_now,1000) == 0 && ~vw(1)
      disp(t_now/1000)
   end

   if n_stim <= 5
      V_all(its,:) = v' ;
      Cai_all(its,:) = cai' ;
   end

   % update CaMK
   CaMKb = CaMKo .* (1.0 - CaMKt) ./ (1.0 + KmCaM ./ cass) ;
   CaMKa = CaMKb + CaMKt ;
   dCaMKt = aCaMK .* CaMKb .* (CaMKb + CaMKt) - bCaMK .* CaMKt ; % Euler step at end

   %% reversal potentials
   ENa = (R .* T ./ F) .* log(nao ./ nai) ;
   EK = (R .* T ./ F) .* log(ko ./ ki) ;
   PKNa = 0.01833 ;
   EKs = (R .* T ./ F) .* log((ko + PKNa .* nao) ./ (ki + PKNa .* nai)) ;

   % convenient shorthand calculations
   vffrt = v .* F .* F ./ (R .* T);
   vfrt = v .* F ./ (R .* T) ;
   frt = F ./ (R .* T) ;

   fINap = (1.0 ./ (1.0 + KmCaMK ./ CaMKa)) ;
   fINaLp = (1.0 ./ (1.0 + KmCaMK ./ CaMKa)) ;
   fItop = (1.0 ./ (1.0 + KmCaMK ./ CaMKa)) ;
   fICaLp = (1.0 ./ (1.0 + KmCaMK ./ CaMKa)) ;

   %% INa (#1)
   [INa, m_, h_, hp_, j_, jp_]...
      = getINa_Grandi(v, m, h, hp, j, jp, fINap, ENa, conductances(:,1), iso, dt, tau_drugs) ;

   %% INaL (#2)
   [INaL, mL_, hL_, hLp_]...
      = getINaL_ORd2011(v, mL, hL, hLp, fINaLp, ENa, celltype, conductances(:,2), dt) ;

   %% Ito (#3)
   [Ito, a_, iF_, iS_, ap_, iFp_, iSp_]...
      = getITo_ORd2011(v, a, iF, iS, ap, iFp, iSp, fItop, EK, celltype, conductances(:,3), dt) ;

   %% ICaL (#4)
   [ICaL_ss, ICaNa_ss, ICaK_ss, ICaL_i, ICaNa_i, ICaK_i, d_, ff_, fs_, fcaf_, fcas_, ...
      jca_, nca_, nca_i_, ffp_, fcafp_, PhiCaL_ss, PhiCaL_i, gammaCaoMyo, gammaCaiMyo]...
      = getICaL_ORd2011_jt(v, d, ff, fs, fcaf, fcas, jca, nca, nca_i, ffp, fcafp, fICaLp,...
      cai, cass, cao, nai, nass, nao, ki, kss, ko, cli, clo, celltype, ICaL_fractionSS, conductances(:,4), iso, dt) ;

   ICaL = ICaL_ss + ICaL_i ;
   ICaNa = ICaNa_ss + ICaNa_i ;
   ICaK = ICaK_ss + ICaK_i ;
   ICaL_tot = ICaL + ICaNa + ICaK ;

   %% IKr (#5)
   [IKr, ikr_c0_, ikr_c1_, ikr_c2_, ikr_o_, ikr_i_] = ...
      getIKr_ORd2011_MM(v, ikr_c0, ikr_c1, ikr_c2, ikr_o, ikr_i, ko, EK, celltype, conductances(:,5), dt) ;

   %% IKs (#6)
   [IKs,xs1_, xs2_]...
      = getIKs_ORd2011(v, xs1, xs2, cai, EKs, celltype, conductances(:,6), iso, dt);

   %% IK1 (#7)
   IK1 = getIK1_CRLP(v, ko, EK, celltype, conductances(:,7));

   %% INaCa (#8)
   [INaCa_i, INaCa_ss]= getINaCa_ORd2011(v, F, R, T, nass, nai, nao, cass, cai, cao, celltype, conductances(:,8), INaCa_fractionSS);

   %% INaK (#9)
   INaK = getINaK_ORd2011(v, F, R, T, nai, nao, ki, ko, celltype, conductances(:,9), iso);

   %% Minor/background currents (#10-13)
   % calculate IKb (#10)
   xkb = 1.0 ./ (1.0 + exp(-(v - 10.8968) ./ 23.9871)) ;

   if ~iso
      GKb = 0.0189 .* conductances(:,10) ;
   else
      GKb = (1 + iso .* (2.5 - 1)) .* 0.0189 .* conductances(:,10) ;
   end

   if celltype == 1
      GKb = GKb .* 0.6 ;
   end
   IKb = GKb .* xkb .* (v - EK) ;

   % calculate INab (#11)
   %PNab = 1.9239e-09 .* conductances(:,11) ;
   INab = PNab .* vffrt .* (nai .* exp(vfrt) - nao) ./ (exp(vfrt) - 1.0) ;

   % calculate ICab (#12)
   % PCab = 5.9194e-08 .* conductances(:,12) ;
   ICab = PCab .* 4.0 .* vffrt .* (gammaCaiMyo .* cai .* exp(2.0 .* vfrt) - gammaCaoMyo .* cao) ./ (exp(2.0 .* vfrt) - 1.0) ;

   % calculate IpCa (#13)
   % GpCa = 5e-04 .* conductances(:,13) ;
   IpCa = GpCa .* cai ./ (0.0005 + cai) ;

   clear xkb GKb % clear temporary variables

   %% Chloride (#14-15)
   % I_ClCa: Ca-activated Cl Current (#14)
   % I_Clbk: background Cl Current (#15)

   % ecl = (R .* T ./ F) .* log(cli ./ clo) ;            % [mV]
   %
   % Fjunc = 1 ;% fraction in SS and in myoplasm - as per literature, I(Ca)Cl is in junctional subspace
   % Fsl = 1 - Fjunc ;
   %
   % GClCa = conductances(:,14) .* 0.2843 ;   % [mS./uF]
   % GClB = conductances(:,15) .* 1.98e-3 ;        % [mS./uF]
   % KdClCa = 0.1 ;    % [mM]

   I_ClCa_junc = Fjunc .* GClCa ./ (1 + KdClCa ./ cass) .* (v - ecl) ;
   I_ClCa_sl = Fsl .* GClCa ./ (1 + KdClCa ./ cai) .* (v - ecl) ;

   I_ClCa_ = I_ClCa_junc + I_ClCa_sl ;
   I_Clbk = GClB .* (v - ecl) ;

   %% Calcium handling - Jrel (#16)
   % calculate ryanodione receptor calcium induced calcium release from the jsr
   fJrelp = (1.0 ./ (1.0 + KmCaMK ./ CaMKa)) ;

   [Jrel, Jrel_np_, Jrel_p_]...
      = getJrel_ORd2011(Jrel_np, Jrel_p, ICaL_ss, cass, cajsr, fJrelp, celltype, conductances(:,16), iso, dt) ;

   %% Jup (#17)
   %calculate serca pump, ca uptake flux
   fJupp = (1.0 ./ (1.0 + KmCaMK ./ CaMKa)) ;
   [Jup, Jleak]...
      = getJup_ORd2011(cai, cansr, fJupp, celltype, conductances(:,17), iso) ;

   % calculate tranlocation flux Jtr
   Jtr = (cansr - cajsr) ./ 60 ;

   %% Iion

   Iion = INa + INaL + Ito + ICaL + ICaNa + ICaK + IKr + IKs + IK1 + INaCa_i...
      + INaCa_ss + INaK + INab + IKb + IpCa + ICab + I_ClCa_ + I_Clbk;

   if t_now < tend
      Istim(stim_cell) = stim_amp * stim_indices(stim_times(:,1) <= t_now & stim_times(:,2) > t_now);
   else
      Istim(stim_cell) = 0;
   end

   Iion = Iion - Istim;

   rhs = Cm/dt*v - Iion ;

   if dt<dt_small*1.0001
      v = A_inv_slow*rhs;
      % A\B is the solution to the equation Ax = B.
      % Matrices A and B must have the same number of rows
      % Same as inv(A)*B, but more accurate
   else
      v = A_inv_fast*rhs;
   end

   % check for stability issues in the simulation and abandon it if
   % necessary
   if isnan(sum(v)) || ~isreal(v)
      warning('stability failure - NaN')
      if ~sum(state_all_last_beat)
         state_all_last_beat...
            = [v, nai, nass, ki, kss, cai, cass, cansr, cajsr, m, hp, h, j,...
            jp, mL, hL, hLp, a, iF, iS, ap, iFp, iSp, d, ff, fs, fcaf, fcas, jca,...
            nca, nca_i, ffp, fcafp, xs1, xs2, Jrel_np, CaMKt, ikr_c0, ikr_c1,...
            ikr_c2, ikr_o, ikr_i, Jrel_p];
         ij_data...
            = [INa, INaL, Ito, ICaL, ICaNa, ICaK, IKr, IKs, IK1, INaCa_i,...
            INaCa_ss, INaK, IKb, INab, ICab, IpCa, I_ClCa_, I_Clbk, Jrel, Jup];
      end
      t = t_now;
      return
   end
   %% calculate diffusion fluxes
   JdiffNa = (nass - nai) ./ 2.0 ;
   JdiffK = (kss - ki) ./ 2.0 ;
   Jdiff = (cass - cai) ./ 0.2 ;

   %% calcium buffer constants

   if ~iso
      kmtrpn = 0.0005 ;
   else
      kmtrpn = (1 + iso .* (1.6 - 1)) .* 0.0005 ;
   end

   %% update intracellular concentrations, using buffers for cai, cass, cajsr
   dnai = -(ICaNa_i + INa + INaL + 3.0 .* INaCa_i + 3.0 .* INaK + INab) .* Acap ./ (F .* vmyo) + JdiffNa .* vss ./ vmyo ;
   dnass = -(ICaNa_ss + 3.0 .* INaCa_ss) .* Acap ./ (F .* vss) - JdiffNa ;

   dki = -(ICaK_i + Ito + IKr + IKs + IK1 + IKb - Istim - 2.0 .* INaK) .* Acap ./ (F .* vmyo) + JdiffK .* vss ./ vmyo ;
   dkss = -ICaK_ss .* Acap ./ (F .* vss) - JdiffK ;

   Bcai = 1.0 ./ (1.0 + cmdnmax .* kmcmdn ./ (kmcmdn + cai) .^ 2.0 + trpnmax .* kmtrpn ./ (kmtrpn + cai) .^ 2.0) ;
   dcai = Bcai .* (-(ICaL_i + IpCa + ICab - 2.0 .* INaCa_i) .* Acap ./ (2.0 .* F .* vmyo) - Jup .* vnsr ./ vmyo + Jdiff .* vss ./ vmyo) ;

   Bcass = 1.0 ./ (1.0 + BSRmax .* KmBSR ./ (KmBSR + cass) .^ 2.0 + BSLmax .* KmBSL ./ (KmBSL + cass) .^ 2.0) ;
   dcass = Bcass .* (-(ICaL_ss - 2.0 .* INaCa_ss) .* Acap ./ (2.0 .* F .* vss) + Jrel .* vjsr ./ vss - Jdiff) ;

   dcansr = Jup - Jtr .* vjsr ./ vnsr ;

   Bcajsr = 1.0 ./ (1.0 + csqnmax .* kmcsqn ./ (kmcsqn + cajsr) .^ 2.0) ;
   dcajsr = Bcajsr .* (Jtr - Jrel) ;

   % EULER
   nai = nai + dnai*dt;
   nass = nass + dnass*dt;
   ki = ki + dki*dt;
   kss = kss + dkss*dt;
   cai = cai + dcai*dt;
   cass = cass + dcass*dt;
   cansr = cansr + dcansr*dt;
   cajsr = cajsr + dcajsr*dt;
   CaMKt = CaMKt + dCaMKt*dt;

   %% update statevars
   m     = m_ ;
   hp    = hp_ ;
   h     = h_ ;
   j     = j_ ;
   jp    = jp_;
   mL    = mL_ ;
   hL    = hL_ ;
   hLp   = hLp_ ;
   a     = a_ ;
   iF    = iF_ ;
   iS    = iS_ ;
   ap    = ap_ ;
   iFp   = iFp_ ;
   iSp   = iSp_ ;
   d     = d_ ;
   ff    = ff_ ;
   fs    = fs_ ;
   fcaf  = fcaf_ ;
   fcas  = fcas_ ;
   jca   = jca_ ;
   nca   = nca_ ;
   nca_i = nca_i_ ;
   ffp   = ffp_ ;
   fcafp = fcafp_ ;
   xs1   = xs1_ ;
   xs2   = xs2_ ;
   Jrel_np = Jrel_np_ ;
   ikr_c0 = ikr_c0_ ;
   ikr_c1 = ikr_c1_ ;
   ikr_c2 = ikr_c2_ ;
   ikr_o = ikr_o_ ;
   ikr_i = ikr_i_ ;
   Jrel_p = Jrel_p_ ;

   if its >= last_beat_startdex && mod(t_now, protocol.downsample_step) == 0
      state_all_last_beat(:, :, round((t_now+protocol.downsample_step)/protocol.downsample_step)-stim_times(end-1,1 ))...
         = [v, nai, nass, ki, kss, cai, cass, cansr, cajsr, m, hp, h, j,...
         jp, mL, hL, hLp, a, iF, iS, ap, iFp, iSp, d, ff, fs, fcaf, fcas, jca,...
         nca, nca_i, ffp, fcafp, xs1, xs2, Jrel_np, CaMKt, ikr_c0, ikr_c1,...
         ikr_c2, ikr_o, ikr_i, Jrel_p] ;
      ij_data(:, :, floor(t_now-t(last_beat_startdex))+1)...
         = [INa, INaL, Ito, ICaL, ICaNa, ICaK, IKr, IKs, IK1, INaCa_i,...
         INaCa_ss, INaK, IKb, INab, ICab, IpCa, I_ClCa_, I_Clbk, Jrel, Jup];
   end

end
%t=t_final;

end
%% LOCAL FUNCTIONS FOR CURRENTS:

%% INa
function [INa, m, h, hp, j, jp] = getINa_Grandi(v, m, h, hp, j, jp, fINap, ENa, INa_Multiplier, iso, dt, tau_drugs)
% The Grandi implementation updated with INa phosphorylation.
%% m gate
mss = 1 ./ ((1 + exp( -(56.86 + v) ./ 9.03 )) .^ 2) ;
taum = 0.1292 .* exp(-((v + 45.79) ./ 15.54) .^ 2) + 0.06487 .* exp(-((v - 4.823) ./ 51.12) .^ 2) ;
%dm = (mss - m) ./ taum ;
taum = tau_drugs.m .* taum;
m = mss-(mss-m)./exp(dt./taum);

%% h gate
ah = (v >= -40) .* (0)...
   + (v < -40) .* (0.057 .* exp( -(v + 80) ./ 6.8)) ;
bh = (v >= -40) .* (0.77 ./ (0.13 .* (1 + exp( -(v + 10.66) ./ 11.1))))...
   + (v < -40) .* ((2.7 .* exp( 0.079 .* v) + 3.1 .* 10 .^ 5 .* exp(0.3485 .* v))) ;
tauh = 1 ./ (ah + bh) ;
tauh = tau_drugs.h .* tauh;

if ~iso
   hss = 1 ./ ((1 + exp( (v + 71.55) ./ 7.43 )) .^ 2) ;
else
   hss = 1 ./ ((1 + exp( (v + 71.55 + (iso .* 5.0)) ./ 7.43 )) .^ 2) ;
end
%dh = (hss - h) ./ tauh ;
h = hss-(hss-h)./exp(dt./tauh); % Rush-Larsen 

%% j gate
aj = (v >= -40) .* (0)...
   + (v < -40) .* (((-2.5428 .* 10 .^ 4 .* exp(0.2444 .* v) - 6.948 .* 10 .^ -6 .* exp(-0.04391 .* v)) .* (v + 37.78)) ./...
   (1 + exp(0.311 .* (v + 79.23)))) ;
bj = (v >= -40) .* ((0.6 .* exp(0.057 .* v)) ./ (1 + exp(-0.1 .* (v + 32))))...
   + (v < -40) .* ((0.02424 .* exp(-0.01052 .* v)) ./ (1 + exp(-0.1378 .* (v + 40.14)))) ;
tauj = 1 ./ (aj + bj) ;
tauj = tau_drugs.j .* tauj;
jss = 1 ./ ((1 + exp((v + 71.55) ./ 7.43)) .^ 2) ;
%dj = (jss - j) ./ tauj ;
j = jss-(jss-j)./exp(dt./tauj);

%% h phosphorylated
if ~iso
   hssp = 1 ./ ((1 + exp((v + 71.55 + 6) ./ 7.43)) .^ 2) ;
else
   hssp = 1 ./ ((1 + exp((v + 71.55 + 6 + (iso .* 5)) ./ 7.43)) .^ 2) ;
end
%dhp = (hssp - hp) ./ tauh ;
hp = hssp-(hssp-hp)./exp(dt./tauh);

%% j phosphorylated
taujp = 1.46 .* tauj ;
%djp = (jss - jp) ./ taujp ;
jp = jss-(jss-jp)./exp(dt./taujp);

if ~iso
   GNa = 11.7802 ;
else
   GNa = (1 + iso .* (2.7 - 1)) .* 11.7802 ;
end
%INa = INa_Multiplier .* GNa .* (v - ENa) .* m .^ 3.0 .* ((1.0 - fINap) .* h .* j + fINap .* hp .* jp) ;
% try make m gate go directly to the steady state
INa = INa_Multiplier .* GNa .* (v - ENa) .* mss .^ 3.0 .* ((1.0 - fINap) .* h .* j + fINap .* hp .* jp) ;
end

%% INaL
function [INaL, mL, hL, hLp] = getINaL_ORd2011(v, mL, hL, hLp, fINaLp, ENa, celltype, INaL_Multiplier, dt)
% calculate INaL
mLss = 1.0 ./ (1.0 + exp((-(v + 42.85)) ./ 5.264)) ;
tm = 0.1292 .* exp(-((v + 45.79) ./ 15.54) .^ 2) + 0.06487 .* exp(-((v - 4.823) ./ 51.12) .^ 2) ;
tmL = tm ;
%dmL = (mLss - mL) ./ tmL ;
mL = mLss-(mLss-mL)./exp(dt./tmL);

hLss = 1.0 ./ (1.0 + exp((v + 87.61) ./ 7.488)) ;
thL = 200.0 ;
%dhL = (hLss - hL) ./ thL ;
hL = hLss-(hLss-hL)./exp(dt./thL);

hLssp = 1.0 ./ (1.0 + exp((v + 93.81) ./ 7.488)) ;
thLp = 3.0 .* thL ;
%dhLp = (hLssp - hLp) ./ thLp ;
hLp = hLssp-(hLssp-hLp)./exp(dt./thLp);

GNaL = 0.0279 .* INaL_Multiplier ;
if celltype == 1
   GNaL = GNaL .* 0.6 ;
end
INaL = GNaL .* (v - ENa) .* mL .* ((1.0 - fINaLp) .* hL + fINaLp .* hLp) ;
end

%% ITo
function [Ito, a, iF, iS, ap, iFp, iSp] = getITo_ORd2011(v, a, iF, iS, ap, iFp, iSp, fItop, EK, celltype, Ito_Multiplier, dt)

% calculate Ito
ass = 1.0 ./ (1.0 + exp((-(v - 14.34)) ./ 14.82)) ;
ta = 1.0515 ./ (1.0 ./ (1.2089 .* (1.0 + exp(-(v - 18.4099) ./ 29.3814))) + 3.5 ./ (1.0 + exp((v + 100.0) ./ 29.3814))) ;
%da = (ass - a) ./ ta ;
a = ass-(ass-a)./exp(dt./ta);

iss = 1.0 ./ (1.0 + exp((v + 43.94) ./ 5.711)) ;
if celltype == 1
   delta_epi = 1.0 - (0.95 ./ (1.0 + exp((v + 70.0) ./ 5.0))) ;
else
   delta_epi = 1.0 ;
end
tiF = 4.562 + 1 ./ (0.3933 .* exp((-(v + 100.0)) ./ 100.0) + 0.08004 .* exp((v + 50.0) ./ 16.59)) ;
tiS = 23.62 + 1 ./ (0.001416 .* exp((-(v + 96.52)) ./ 59.05) + 1.780e-8 .* exp((v + 114.1) ./ 8.079)) ;
tiF = tiF .* delta_epi ;
tiS = tiS .* delta_epi ;
AiF = 1.0 ./ (1.0 + exp((v - 213.6) ./ 151.2)) ;
AiS = 1.0 - AiF ;
%diF = (iss - iF) ./ tiF ;
%diS = (iss - iS) ./ tiS ;
iF = iss-(iss-iF)./exp(dt./tiF);
iS = iss-(iss-iS)./exp(dt./tiS);
i = AiF .* iF + AiS .* iS ;

assp = 1.0 ./ (1.0 + exp((-(v - 24.34)) ./ 14.82)) ;
%dap = (assp - ap) ./ ta ;
ap = assp-(assp-ap)./exp(dt./ta);

dti_develop = 1.354 + 1.0e-4 ./ (exp((v - 167.4) ./ 15.89) + exp(-(v - 12.23) ./ 0.2154)) ;
dti_recover = 1.0 - 0.5 ./ (1.0 + exp((v + 70.0) ./ 20.0)) ;
tiFp = dti_develop .* dti_recover .* tiF ;
tiSp = dti_develop .* dti_recover .* tiS ;
% diFp = (iss - iFp) ./ tiFp ;
% diSp = (iss - iSp) ./ tiSp ;
iFp = iss-(iss-iFp)./exp(dt./tiFp);
iSp = iss-(iss-iSp)./exp(dt./tiSp);
ip = AiF .* iFp + AiS .* iSp ;

Gto = 0.16 .* Ito_Multiplier ;
if celltype == 1
   Gto = Gto .* 2.0 ;
elseif celltype == 2
   Gto = Gto .* 2.0 ;
end

Ito = Gto .* (v - EK) .* ((1.0 - fItop) .* a .* i + fItop .* ap .* ip) ;
end

%% ICaL
% a variant updated by J. Tomek, using a changed activation curve
% it computes both ICaL in subspace and myoplasm (_i)
function [ICaL_ss, ICaNa_ss, ICaK_ss, ICaL_i, ICaNa_i, ICaK_i, d, ff, fs,...
   fcaf, fcas, jca, nca, nca_i, ffp, fcafp, PhiCaL_ss, PhiCaL_i,...
   gammaCaoMyo, gammaCaiMyo]...
   = getICaL_ORd2011_jt(v, d, ff, fs, fcaf, fcas, jca, nca, nca_i, ffp, fcafp,...
   fICaLp, cai, cass, cao, nai, nass, nao, ki, kss, ko, cli, clo, celltype, ICaL_fractionSS, ICaL_PCaMultiplier, iso, dt)

% physical constants
R = 8314.0 ;
T = 310.0 ;
F = 96485.0 ;
vffrt = v .* F .* F ./ (R .* T) ;
vfrt = v .* F ./ (R .* T) ;

% calculate ICaL, ICaNa, ICaK

if ~iso
   dss = 1.0763 .* exp(-1.0070 .* exp(-0.0829 .* v)) ;  % magyar
else
   dss = 1.0763 .* exp(-1.0070 .* exp(-0.0829 .* (v + (iso .* 12.0)))) ;  % magyar
end
if(v > 31.4978) % activation cannot be greater than 1
   dss = 1;
end
td = 0.6 + 1.0 ./ (exp(-0.05 .* (v + 6.0)) + exp(0.09 .* (v + 14.0))) ;
%dd = (dss - d) ./ td ;
d = dss-(dss-d)./exp(dt./td);

if ~iso
   fss = 1.0 ./ (1.0 + exp((v + 19.58) ./ 3.696)) ;
else
   fss = 1.0 ./ (1.0 + exp((v + 19.58 + (iso .* 8.0)) ./ 3.696)) ;
end

tff = 7.0 + 1.0 ./ (0.0045 .* exp(-(v + 20.0) ./ 10.0) + 0.0045 .* exp((v + 20.0) ./ 10.0)) ;
tfs = 1000.0 + 1.0 ./ (0.000035 .* exp(-(v + 5.0) ./ 4.0) + 0.000035 .* exp((v + 5.0) ./ 6.0)) ;
Aff = 0.6 ;
Afs = 1.0 - Aff ;
% dff = (fss - ff) ./ tff ;
% dfs = (fss - fs) ./ tfs ;
ff = fss-(fss-ff)./exp(dt./tff);
fs = fss-(fss-fs)./exp(dt./tfs);
f = Aff .* ff + Afs .* fs ;

fcass = fss ;
tfcaf = 7.0 + 1.0 ./ (0.04 .* exp(-(v - 4.0) ./ 7.0) + 0.04 .* exp((v - 4.0) ./ 7.0)) ;
tfcas = 100.0 + 1.0 ./ (0.00012 .* exp(-v ./ 3.0) + 0.00012 .* exp(v ./ 7.0)) ;
Afcaf = 0.3 + 0.6 ./ (1.0 + exp((v - 10.0) ./ 10.0)) ;
Afcas = 1.0 - Afcaf ;
% dfcaf = (fcass - fcaf) ./ tfcaf ;
% dfcas = (fcass - fcas) ./ tfcas ;
fcaf = fcass-(fcass-fcaf)./exp(dt./tfcaf);
fcas = fcass-(fcass-fcas)./exp(dt./tfcas);
fca = Afcaf .* fcaf + Afcas .* fcas ;

tjca = 75 ;
jcass = 1.0 ./ (1.0 + exp((v + 18.08) ./ 2.7916)) ;
%djca = (jcass - jca) ./ tjca ;
jca = jcass-(jcass-jca)./exp(dt./tjca);

tffp = 2.5 .* tff ;
%dffp = (fss - ffp) ./ tffp ;
ffp = fss-(fss-ffp)./exp(dt./tffp);

fp = Aff .* ffp + Afs .* fs ;
tfcafp = 2.5 .* tfcaf ;
%dfcafp = (fcass - fcafp) ./ tfcafp ;
fcafp = fcass-(fcass-fcafp)./exp(dt./tfcafp);
fcap = Afcaf .* fcafp + Afcas .* fcas ;

% SS nca
Kmn = 0.002 ;
k2n = 500.0 ;
km2n = jca .* 1 ;
anca = 1.0 ./ (k2n ./ km2n + (1.0 + Kmn ./ cass) .^ 4.0) ;
dnca = anca .* k2n - nca .* km2n ;
% EULER
nca = nca + dnca .* dt;

% myoplasmic nca
anca_i = 1.0 ./ (k2n ./ km2n + (1.0 + Kmn ./ cai) .^ 4.0) ;
dnca_i = anca_i .* k2n - nca_i .* km2n ;
% EULER
nca_i = nca_i + dnca_i .* dt;

% SS driving force
Io = 0.5 .* (nao + ko + clo + 4 .* cao) ./ 1000 ; % ionic strength outside. ./1000 is for things being in micromolar
Ii = 0.5 .* (nass + kss + cli + 4 .* cass) ./ 1000 ; % ionic strength outside. ./1000 is for things being in micromolar
% The ionic strength is too high for basic DebHuc. We'll use Davies
dielConstant = 74 ; % water at 37°.
temp = 310 ; % body temp in kelvins.
constA = 1.82 .* 10 .^ 6 .* (dielConstant .* temp) .^ (-1.5) ;

gamma_cai = exp(-constA .* 4 .* (sqrt(Ii) ./ (1 + sqrt(Ii)) - 0.3 .* Ii)) ;
gamma_cao = exp(-constA .* 4 .* (sqrt(Io) ./ (1 + sqrt(Io)) - 0.3 .* Io)) ;
gamma_nai = exp(-constA .* 1 .* (sqrt(Ii) ./ (1 + sqrt(Ii)) - 0.3 .* Ii)) ;
gamma_nao = exp(-constA .* 1 .* (sqrt(Io) ./ (1 + sqrt(Io)) - 0.3 .* Io)) ;
gamma_ki = exp(-constA .* 1 .* (sqrt(Ii) ./ (1 + sqrt(Ii)) - 0.3 .* Ii)) ;
gamma_kao = exp(-constA .* 1 .* (sqrt(Io) ./ (1 + sqrt(Io))-  0.3 .* Io)) ;


PhiCaL_ss = 4.0 .* vffrt .* (gamma_cai .* cass .* exp(2.0 .* vfrt) - gamma_cao .* cao) ./ (exp(2.0 .* vfrt) - 1.0) ;
PhiCaNa_ss = 1.0 .* vffrt .* (gamma_nai .* nass .* exp(1.0 .* vfrt) - gamma_nao .* nao) ./ (exp(1.0 .* vfrt) - 1.0) ;
PhiCaK_ss = 1.0 .* vffrt .* (gamma_ki .* kss .* exp(1.0 .* vfrt) - gamma_kao .* ko) ./ (exp(1.0 .* vfrt) - 1.0) ;

% Myo driving force
Io = 0.5 .* (nao + ko + clo + 4 .* cao) ./ 1000 ; % ionic strength outside. ./1000 is for things being in micromolar
Ii = 0.5 .* (nai + ki + cli + 4 .* cai) ./ 1000 ; % ionic strength outside. ./1000 is for things being in micromolar
% The ionic strength is too high for basic DebHuc. We'll use Davies
dielConstant = 74 ; % water at 37°.
temp = 310 ; % body temp in kelvins.
constA = 1.82 .* 10 .^ 6 .* (dielConstant .* temp) .^ (-1.5) ;

gamma_cai = exp(-constA .* 4 .* (sqrt(Ii) ./ (1 + sqrt(Ii)) - 0.3 .* Ii)) ;
gamma_cao = exp(-constA .* 4 .* (sqrt(Io) ./ (1 + sqrt(Io)) - 0.3 .* Io)) ;
gamma_nai = exp(-constA .* 1 .* (sqrt(Ii) ./ (1 + sqrt(Ii)) - 0.3 .* Ii)) ;
gamma_nao = exp(-constA .* 1 .* (sqrt(Io) ./ (1 + sqrt(Io)) - 0.3 .* Io)) ;
gamma_ki = exp(-constA .* 1 .* (sqrt(Ii) ./ (1 + sqrt(Ii)) - 0.3 .* Ii)) ;
gamma_kao = exp(-constA .* 1 .* (sqrt(Io) ./ (1 + sqrt(Io)) - 0.3 .* Io)) ;

gammaCaoMyo = gamma_cao ;
gammaCaiMyo = gamma_cai ;

PhiCaL_i = 4.0 .* vffrt .* (gamma_cai .* cai .* exp(2.0 .* vfrt) - gamma_cao .* cao) ./ (exp(2.0 .* vfrt) - 1.0) ;
PhiCaNa_i = 1.0 .* vffrt .* (gamma_nai .* nai .* exp(1.0 .* vfrt) - gamma_nao .* nao) ./ (exp(1.0 .* vfrt) - 1.0) ;
PhiCaK_i = 1.0 .* vffrt .* (gamma_ki .* ki .* exp(1.0 .* vfrt) - gamma_kao .* ko) ./ (exp(1.0 .* vfrt) - 1.0) ;

% The rest
if ~iso
   PCa = 8.3757e-05 .* ICaL_PCaMultiplier ;
else
   PCa = 8.3757e-05 .* ICaL_PCaMultiplier .* (1 + iso .* (1.9 - 1)) ;
end

if celltype == 1
   PCa = PCa .* 1.2 ;
elseif celltype == 2
   PCa = PCa .* 2 ;
end

PCap = 1.1 .* PCa ;
PCaNa = 0.00125 .* PCa ;
PCaK = 3.574e-4 .* PCa ;
PCaNap = 0.00125 .* PCap ;
PCaKp = 3.574e-4 .* PCap ;

ICaL_ss = (1.0 - fICaLp) .* PCa .* PhiCaL_ss .* d .* (f .* (1.0 - nca) + jca .* fca .* nca)...
   + fICaLp .* PCap .* PhiCaL_ss.*  d .* (fp .* (1.0 - nca) + jca .* fcap .* nca) ;
ICaNa_ss = (1.0 - fICaLp) .* PCaNa .* PhiCaNa_ss .* d .* (f .* (1.0 - nca) + jca .* fca .* nca)...
   + fICaLp .* PCaNap .* PhiCaNa_ss .* d .* (fp .* (1.0 - nca) + jca .* fcap .* nca) ;
ICaK_ss = (1.0 - fICaLp) .* PCaK .* PhiCaK_ss .* d .* (f .* (1.0 - nca) + jca .* fca .* nca)...
   + fICaLp .* PCaKp .* PhiCaK_ss .* d .* (fp .* (1.0 - nca) + jca .* fcap .* nca) ;

ICaL_i = (1.0 - fICaLp) .* PCa .* PhiCaL_i .* d .* (f .* (1.0 - nca_i) + jca .* fca .* nca_i)...
   + fICaLp .* PCap .* PhiCaL_i .* d .* (fp .* (1.0 - nca_i) + jca .* fcap .* nca_i) ;
ICaNa_i = (1.0 - fICaLp) .* PCaNa .* PhiCaNa_i .* d .* (f .* (1.0 - nca_i) + jca .* fca .* nca_i)...
   + fICaLp .* PCaNap .* PhiCaNa_i .* d .* (fp .* (1.0 - nca_i) + jca .* fcap .* nca_i) ;
ICaK_i = (1.0 - fICaLp) .* PCaK .* PhiCaK_i .* d .* (f .* (1.0 - nca_i) + jca .* fca .* nca_i)...
   + fICaLp .* PCaKp .* PhiCaK_i .* d .* (fp .* (1.0 - nca_i) + jca .* fcap .* nca_i) ;

% And we weight ICaL (in ss) and ICaL_i
ICaL_i = ICaL_i .* (1 - ICaL_fractionSS) ;
ICaNa_i = ICaNa_i .* (1 - ICaL_fractionSS) ;
ICaK_i = ICaK_i .* (1 - ICaL_fractionSS) ;
ICaL_ss = ICaL_ss .* ICaL_fractionSS ;
ICaNa_ss = ICaNa_ss .* ICaL_fractionSS ;
ICaK_ss = ICaK_ss .* ICaL_fractionSS ;

end

%% IKr
% Variant based on Lu-Vandenberg
function [IKr, c0, c1, c2, o, i ]...
   = getIKr_ORd2011_MM(V, c0, c1, c2, o, i, ko, EK, celltype, IKr_Multiplier, dt)

% physical constants
R = 8314.0 ;
T = 310.0 ;
F = 96485.0 ;

b = 0; % no channels blocked in via the mechanism of specific MM states
vfrt = V .* F ./ (R .* T) ;

% transition rates
% from c0 to c1 in l-v model,
alpha = 0.1161 .* exp(0.2990 .* vfrt) ;
% from c1 to c0 in l-v./
beta = 0.2442 .* exp(-1.604 .* vfrt) ;

% from c1 to c2 in l-v./
alpha1 = 1.25 .* 0.1235 ;
% from c2 to c1 in l-v./
beta1 = 0.1911 ;

% from c2 to o./           c1 to o
alpha2 = 0.0578 .* exp(0.9710 .* vfrt) ;
% from o to c2./
beta2 = 0.349e-3 .* exp(-1.062 .* vfrt) ;

% from o to i
alphai = 0.2533 .* exp(0.5953 .* vfrt) ;
% from i to o
betai = 1.25 .* 0.0522 .* exp(-0.8209 .* vfrt) ;

% from c2 to i (from c1 in orig)
alphac2ToI = 0.52e-4 .* exp(1.525 .* vfrt) ;
% from i to c2
% betaItoC2 = 0.85e-8 .* exp(-1.842 .* vfrt);
betaItoC2 = (beta2 .* betai .* alphac2ToI) ./ (alpha2 .* alphai) ;
% transitions themselves
% for reason of backward compatibility of naming of an older version of a
% MM IKr, c3 in code is c0 in article diagram, c2 is c1, c1 is c2.

dc0 = c1 .* beta - c0 .* alpha ; % delta for c0
dc1 = c0 .* alpha + c2 .* beta1 - c1 .* (beta + alpha1) ; % c1
dc2 = c1 .* alpha1 + o .* beta2 + i .* betaItoC2 - c2 .* (beta1 + alpha2 + alphac2ToI) ; % subtraction is into c2, to o, to i. % c2
do = c2 .* alpha2 + i .* betai - o .* (beta2 + alphai);
di = c2 .* alphac2ToI + o .* alphai - i .* (betaItoC2 + betai) ;

% EULER
c0 = c0 + dc0*dt;
c1 = c1 + dc1*dt;
c2 = c2 + dc2*dt;
o = o + do*dt;
i = i + di*dt;

GKr = 0.0321 .* sqrt(ko ./ 5) .* IKr_Multiplier ; % 1st element compensates for change to ko (sqrt(5./5.4).* 0.0362)
if celltype == 1
   GKr = GKr .* 1.3 ;
elseif celltype == 2
   GKr = GKr .* 0.8 ;
end

IKr = GKr .* o .* (V - EK) ;
end

%% IKs
function [IKs,xs1, xs2]...
   = getIKs_ORd2011(v, xs1, xs2, cai, EKs, celltype, IKs_Multiplier, iso, dt)
% calculate IKs
xs1ss = 1.0 ./ (1.0 + exp((-(v + 11.60)) ./ 8.932)) ;

if ~iso
   txs1 = 817.3 + 1.0 ./ (2.326e-4 .* exp((v + 48.28) ./ 17.80) + 0.001292 .* exp((-(v + 210.0)) ./ 230.0)) ;
else
   txs1 = 817.3...
      - (iso .* 1.75) .* (1.0 ./ (2.326e-4 .* exp((10 + 48.28) ./ 17.80) + 0.001292 .* exp((-(10 + 210.0)) ./ 230.0)))...
      + (1 + iso .* (2.75 - 1)) .* (1.0 ./ (2.326e-4 .* exp((v + 48.28) ./ 17.80) + 0.001292 .* exp((-(v + 210.0)) ./ 230.0)));
end
% dxs1 = (xs1ss - xs1) ./ txs1 ;
xs1 = xs1ss-(xs1ss-xs1)./exp(dt./txs1);

xs2ss = xs1ss ;
txs2 = 1.0 ./ (0.01 .* exp((v - 50.0) ./ 20.0) + 0.0193 .* exp((-(v + 66.54)) ./ 31.0)) ;
% dxs2 = (xs2ss - xs2) ./ txs2 ;
xs2 = xs2ss-(xs2ss-xs2)./exp(dt./txs2);

KsCa = 1.0 + 0.6 ./ (1.0 + (3.8e-5 ./ cai) .^ 1.4) ;

if ~iso
   GKs = 0.0011 .* IKs_Multiplier ;
else
   GKs = (1 + iso .* (8.0 - 1)) .* 0.0011 .* IKs_Multiplier ;
end

if celltype == 1
   GKs = GKs .* 1.4 ;
end
IKs = GKs .* KsCa .* xs1 .* xs2 .* (v - EKs) ;
end

%% IK1
function [IK1]...
   = getIK1_CRLP(v, ko, EK, celltype, IK1_Multiplier)
% IK1
aK1 = 4.094 ./ (1 + exp(0.1217 .* (v - EK - 49.934))) ;
bK1 = (15.72 .* exp(0.0674 .* (v - EK - 3.257)) + exp(0.0618 .* (v - EK - 594.31))) ./ (1 + exp(-0.1629 .* (v - EK + 14.207))) ;
K1ss = aK1 ./ (aK1 + bK1) ;

GK1 = IK1_Multiplier .* 0.6992 ; %0.7266; %.* sqrt(5./5.4))
if celltype == 1
   GK1 = GK1 .* 1.2 ;
elseif celltype == 2
   GK1 = GK1 .* 1.3 ;
end
IK1 = GK1 .* sqrt(ko ./ 5) .* K1ss .* (v - EK) ;
end

%% INaCa
function [INaCa_i, INaCa_ss]...
   = getINaCa_ORd2011(v, F, R, T, nass, nai, nao, cass, cai, cao, celltype, INaCa_Multiplier, INaCa_fractionSS)
zca = 2.0 ;
kna1 = 15.0 ;
kna2 = 5.0 ;
kna3 = 88.12 ;
kasymm = 12.5 ;
wna = 6.0e4 ;
wca = 6.0e4 ;
wnaca = 5.0e3 ;
kcaon = 1.5e6 ;
kcaoff = 5.0e3 ;
qna = 0.5224 ;
qca = 0.1670 ;
hca = exp((qca .* v .* F) ./ (R .* T)) ;
hna = exp((qna .* v .* F) ./ (R .* T)) ;
h1 = 1 + nai ./ kna3 .* (1 + hna) ;
h2 = (nai .* hna) ./ (kna3 .* h1) ;
h3 = 1.0 ./ h1 ;
h4 = 1.0 + nai ./ kna1 .* (1 + nai ./ kna2) ;
h5 = nai .* nai ./ (h4 .* kna1 .* kna2) ;
h6 = 1.0 ./ h4 ;
h7 = 1.0 + nao ./ kna3 .* (1.0 + 1.0 ./ hna) ;
h8 = nao ./ (kna3 .* hna .* h7) ;
h9 = 1.0 ./ h7 ;
h10 = kasymm + 1.0 + nao ./ kna1 .* (1.0 + nao ./ kna2) ;
h11 = nao .* nao ./ (h10 .* kna1 .* kna2) ;
h12 = 1.0 ./ h10 ;
k1 = h12 .* cao .* kcaon ;
k2 = kcaoff ;
k3p = h9 .* wca ;
k3pp = h8 .* wnaca ;
k3 = k3p + k3pp ;
k4p = h3 .* wca ./ hca ;
k4pp = h2 .* wnaca ;
k4 = k4p + k4pp ;
k5 = kcaoff ;
k6 = h6 .* cai .* kcaon ;
k7 = h5 .* h2 .* wna ;
k8 = h8 .* h11 .* wna ;
x1 = k2 .* k4 .* (k7 + k6) + k5 .* k7 .* (k2 + k3) ;
x2 = k1 .* k7 .* (k4 + k5) + k4 .* k6 .* (k1 + k8) ;
x3 = k1 .* k3 .* (k7 + k6) + k8 .* k6 .* (k2 + k3) ;
x4 = k2 .* k8 .* (k4 + k5) + k3 .* k5 .* (k1 + k8) ;
E1 = x1 ./ (x1 + x2 + x3 + x4) ;
E2 = x2 ./ (x1 + x2 + x3 + x4) ;
E3 = x3 ./ (x1 + x2 + x3 + x4) ;
E4 = x4 ./ (x1 + x2 + x3 + x4) ;
KmCaAct = 150.0e-6 ;
allo = 1.0 ./ (1.0 + (KmCaAct ./ cai) .^ 2.0) ;
zna = 1.0 ;
JncxNa = 3.0 .* (E4 .* k7 - E1 .* k8) + E3 .* k4pp - E2 .* k3pp ;
JncxCa = E2 .* k2 - E1 .* k1 ;
Gncx = 0.0034 .* INaCa_Multiplier ;
if celltype == 1
   Gncx = Gncx .* 1.1 ;
elseif celltype == 2
   Gncx = Gncx .* 1.4 ;
end
INaCa_i = (1 - INaCa_fractionSS) .* Gncx .* allo .* (zna .* JncxNa + zca .* JncxCa) ;

% calculate INaCa_ss
h1 = 1 + nass ./ kna3 .* (1 + hna) ;
h2 = (nass .* hna) ./ (kna3 .* h1) ;
h3 = 1.0 ./ h1 ;
h4 = 1.0 + nass ./ kna1 .* (1 + nass ./ kna2) ;
h5 = nass .* nass ./ (h4 .* kna1 .* kna2) ;
h6 = 1.0 ./ h4 ;
h7 = 1.0 + nao ./ kna3 .* (1.0 + 1.0 ./ hna) ;
h8 = nao ./ (kna3 .* hna .* h7) ;
h9 = 1.0 ./ h7 ;
h10 = kasymm + 1.0 + nao ./ kna1 .* (1 + nao ./ kna2) ;
h11 = nao .* nao ./ (h10 .* kna1 .* kna2) ;
h12 = 1.0 ./ h10 ;

k1 = h12 .* cao .* kcaon ;
k2 = kcaoff ;
k3p = h9 .* wca ;
k3pp = h8 .* wnaca ;
k3 = k3p + k3pp ;
k4p = h3 .* wca ./ hca ;
k4pp = h2 .* wnaca ;
k4 = k4p + k4pp ;
k5 = kcaoff ;
k6 = h6 .* cass .* kcaon ;
k7 = h5 .* h2 .* wna ;
k8 = h8 .* h11 .* wna ;
x1 = k2 .* k4 .* (k7 + k6) + k5 .* k7 .* (k2 + k3) ;
x2 = k1 .* k7 .* (k4 + k5) + k4 .* k6 .* (k1 + k8) ;
x3 = k1 .* k3 .* (k7 + k6) + k8 .* k6 .* (k2 + k3) ;
x4 = k2 .* k8 .* (k4 + k5) + k3 .* k5 .* (k1 + k8) ;
E1 = x1 ./ (x1 + x2 + x3 + x4) ;
E2 = x2 ./ (x1 + x2 + x3 + x4) ;
E3 = x3 ./ (x1 + x2 + x3 + x4) ;
E4 = x4 ./ (x1 + x2 + x3 + x4) ;
KmCaAct = 150.0e-6 ;
allo = 1.0 ./ (1.0 + (KmCaAct ./ cass) .^ 2.0) ;
JncxNa = 3.0 .* (E4 .* k7 - E1 .* k8) + E3 .* k4pp - E2 .* k3pp ;
JncxCa = E2 .* k2 - E1 .* k1 ;
INaCa_ss = INaCa_fractionSS .* Gncx .* allo .* (zna .* JncxNa + zca .* JncxCa) ;
end

%% INaK
function INaK = getINaK_ORd2011(v, F, R, T, nai, nao, ki, ko, celltype, INaK_Multiplier, iso)
% calculate INaK
zna = 1.0 ;
k1p = 949.5 ;
k1m = 182.4 ;
k2p = 687.2 ;
k2m = 39.4 ;
k3p = 1899.0 ;
k3m = 79300.0 ;
k4p = 639.0 ;
k4m = 40.0 ;

if ~iso
   Knai0 = 9.073 ;
else
   Knai0 = (1 + iso .* (0.7 - 1)) .* 9.073 ;
end

Knao0 = 27.78 ;
delta = -0.1550 ;
Knai = Knai0 .* exp((delta .* v .* F) ./ (3.0 .* R .* T)) ;
Knao = Knao0 .* exp(((1.0 - delta) .* v .* F) ./ (3.0 .* R .* T)) ;
Kki = 0.5 ;
Kko = 0.3582 ;
MgADP = 0.05 ;
MgATP = 9.8 ;
Kmgatp = 1.698e-7 ;
H = 1.0e-7 ;
eP = 4.2 ;
Khp = 1.698e-7 ;
Knap = 224.0 ;
Kxkur = 292.0 ;
P = eP ./ (1.0 + H ./ Khp + nai ./ Knap + ki ./ Kxkur) ;
a1 = (k1p .* (nai ./ Knai) .^ 3.0) ./ ((1.0 + nai ./ Knai) .^ 3.0 + (1.0 + ki ./ Kki) .^ 2.0 - 1.0) ;
b1 = k1m .* MgADP ;
a2 = k2p ;
b2 = (k2m .* (nao ./ Knao) .^ 3.0) ./ ((1.0 + nao ./ Knao) .^ 3.0 + (1.0 + ko ./ Kko) .^ 2.0 - 1.0) ;
a3 = (k3p .* (ko ./ Kko) .^ 2.0) ./ ((1.0 + nao ./ Knao) .^ 3.0 + (1.0 + ko ./ Kko) .^ 2.0 - 1.0) ;
b3 = (k3m .* P .* H) ./ (1.0 + MgATP ./ Kmgatp) ;
a4 = (k4p .* MgATP ./ Kmgatp) ./ (1.0 + MgATP ./ Kmgatp) ;
b4 = (k4m .* (ki ./ Kki) .^ 2.0) ./ ((1.0 + nai ./ Knai) .^ 3.0 + (1.0 + ki ./ Kki) .^ 2.0 - 1.0) ;
x1 = a4 .* a1 .* a2 + b2 .* b4 .* b3 + a2 .* b4 .* b3 + b3 .* a1 .* a2 ;
x2 = b2 .* b1 .* b4 + a1 .* a2 .* a3 + a3 .* b1 .* b4 + a2 .* a3 .* b4 ;
x3 = a2 .* a3 .* a4 + b3 .* b2 .* b1 + b2 .* b1 .* a4 + a3 .* a4 .* b1 ;
x4 = b4 .* b3 .* b2 + a3 .* a4 .* a1 + b2 .* a4 .* a1 + b3 .* b2 .* a1 ;
E1 = x1 ./ (x1 + x2 + x3 + x4) ;
E2 = x2 ./ (x1 + x2 + x3 + x4) ;
E3 = x3 ./ (x1 + x2 + x3 + x4) ;
E4 = x4 ./ (x1 + x2 + x3 + x4) ;
zk = 1.0 ;
JnakNa = 3.0 .* (E1 .* a3 - E2 .* b3) ;
JnakK = 2.0 .* (E4 .* b1 - E3 .* a1) ;
Pnak = 15.4509 .* INaK_Multiplier ;
if celltype == 1
   Pnak = Pnak .* 0.9 ;
elseif celltype == 2
   Pnak = Pnak .* 0.7 ;
end
INaK = Pnak .* (zna .* JnakNa + zk .* JnakK) ;
end

%% Jrel
function [Jrel, Jrelnp, Jrelp]...
   = getJrel_ORd2011(Jrelnp, Jrelp, ICaL, cass, cajsr, fJrelp, celltype, Jrel_Multiplier, iso, dt)

jsrMidpoint = 1.7 ;
bt = 4.75 ;

if ~iso
   a_rel = 0.5 .* bt ;
else
   a_rel = (1 + iso .* (2.5 - 1)) .* 0.5 .* bt ;
end
Jrel_inf = a_rel .* (-ICaL) ./ (1.0 + (jsrMidpoint ./ cajsr) .^ 8.0) ;
if celltype == 2
   Jrel_inf = Jrel_inf .* 1.7 ;
end

if ~iso
   tau_rel = bt ./ (1.0 + 0.0123 ./ cajsr) ;
else
   tau_rel = (1 + iso .* (0.75 - 1)) .* (bt ./ (1.0 + 0.0123 ./ cajsr)) ;
end

if tau_rel < 0.001
   tau_rel = 0.001 ;
end

%dJrelnp = (Jrel_inf - Jrelnp) ./ tau_rel ;
Jrelnp = Jrel_inf-(Jrel_inf-Jrelnp)./exp(dt./tau_rel);

btp = 1.25 .* bt ;

if ~iso
   a_relp = 0.5 .* btp ;
else
   a_relp = (1 + iso .* (2.5 - 1)) .* 0.5 .* btp ;
end

Jrel_infp = a_relp .* (-ICaL) ./ (1.0 + (jsrMidpoint ./ cajsr) .^ 8.0) ;
if celltype == 2
   Jrel_infp = Jrel_infp .* 1.7 ;
end

if ~iso
   tau_relp = btp ./ (1.0 + 0.0123 ./ cajsr) ;
else
   tau_relp = (1 + iso .* (0.75 - 1)) .* (btp ./ (1.0 + 0.0123 ./ cajsr)) ;
end

if tau_relp < 0.001
   tau_relp = 0.001 ;
end

%dJrelp = (Jrel_infp - Jrelp) ./ tau_relp ;
Jrelp = Jrel_infp-(Jrel_infp-Jrelp)./exp(dt./tau_relp);

Jrel = Jrel_Multiplier .* 1.5378 .* ((1.0 - fJrelp) .* Jrelnp + fJrelp .* Jrelp) ;
end

%% Jup
function [Jup, Jleak]...
   = getJup_ORd2011(cai, cansr, fJupp, celltype, Jup_Multiplier, iso)

%calculate serca pump, ca uptake flux
if ~iso
   Jupnp = Jup_Multiplier .* 0.005425 .* cai ./ (cai + 0.00092) ;
else
   Jupnp = Jup_Multiplier .* 0.005425 .* cai ./ (cai + (1 + iso .* (0.54 - 1)) .* 0.00092) ;
end

if ~iso
   Jupp = Jup_Multiplier .* 2.75 .* 0.005425 .* cai ./ (cai + 0.00092 - 0.00017) ;
else
   Jupp = Jup_Multiplier .* 2.75 .* 0.005425 .* cai ./ (cai + (1 + iso .* (0.54 - 1)) .* (0.00092 - 0.00017)) ;
end

if celltype == 1
   Jupnp = Jupnp .* 1.3 ;
   Jupp = Jupp .* 1.3 ;
end

Jleak = Jup_Multiplier .* 0.0048825 .* cansr ./ 15.0 ;
Jup = (1.0 - fJupp) .* Jupnp + fJupp .* Jupp - Jleak ;
end





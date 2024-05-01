function [] = VW(vw, tissue, sample, minerva, drugs, ascend, rj_s2a)
% The first argument vw is a 4-element vector as follows:
% vw(1) = flag that this simulation will be calculating VW,
% vw(2) = S2 position as a fraction of cable length
% vw(3,4) = initial start and end points for S1-S2 intervals (starting "window")

%% Intialize Setup
% initialize variables that may have been excluded
if ~exist('minerva', 'var')
   minerva = 0;
end

if ~exist('drugs', 'var')
   drugs.effect = 1;
   drugs.target = 1;
end

N = 100 ; % 100 CMs in cable
protocol = struct(); % will hold details of protocol
s2cellnum = floor(vw(2) * N * 1.00001); % S2-CM position
vw_failure = 0; % flag if VW could not be calculated

rj = rj_s2a(1); % Rj multiplier
s2a = rj_s2a(2); % S2 amplitude

protocol.iso = 0; % toggle to turn on isoproterenol application, 0=off, 1=on

protocol.Rj_mod = rj;
% generate a tag for dynamic filename creation
if rj ~=5
   protocol.rjtag = ['_Rj_', num2str(rj)];
else
   protocol.rjtag = '';
end

protocol.ascend = ascend;
% generate a tag for dynamic filename creation
if protocol.ascend > 0
   protocol.ascendtag = ['_ascendingG', num2str(protocol.ascend)];
elseif protocol.ascend < 0
   protocol.ascendtag = ['_descendingG', num2str(protocol.ascend*-1)];
else
   protocol.ascendtag = '';
end

protocol.drug_effect = drugs.effect; % eg. 0.7, 1.3
protocol.drug_target = drugs.target; % # of conductance
protocol.drug_loc = drugs.drug_loc; % location in cable of drug

protocol.drugtag = ''; % generate a tag for dynamic filename creation
if protocol.drug_effect ~= 1
   if  drugs.drug_loc(1) == 1 && drugs.drug_loc(2) == 100
   drug_loc = '';
   elseif drugs.drug_loc(2) == 49 && drugs.drug_loc(1) == 1
      drug_loc = '_upstream';
   elseif drugs.drug_loc(1) == 51 && drugs.drug_loc(2) == 100
      drug_loc = '_downstream';
   else 
      drug_loc = ['_CM',num2str(drugs.drug_loc(1)),'-',num2str(drugs.drug_loc(2))];
   end
   protocol.drugtag = ['_G',num2str(protocol.drug_target),'_',num2str(protocol.drug_effect), drug_loc];
end

% uniform stretch of cell #50
protocol.stretch = drugs.unistretch;
% generate a tag for dynamic filename creation
if length(drugs.unistretch) > 1
protocol.stretchtag = ['_stretch', num2str(length(protocol.stretch))];
else
   protocol.stretchtag = '';
end

% make tissue numbering strings for creating filenames
if tissue > 1
   tissue_string = num2str(tissue);
else
   % by default all uniform cables are created from a cell drawn from a 
   % population with sigma = 0.1. This can be changed by setting this 
   % variable, and its partner in f_tomek_cable_CNRL_CVTS.m (line 80) to
   % a different value
   uniform_variant = 3; 
   tissue_string = num2str(tissue + uniform_variant * 0.1);
end

%% Establish initial conditions (ICs) if not already done

% create filenames
ic_filename_partial = ['Rush-Larsen/ICs/ics_N', num2str(N),...
   '_t', tissue_string,...
   '_s', num2str(sample)];
ic_filename=[ic_filename_partial, '_beat21', protocol.drugtag, protocol.ascendtag, protocol.rjtag, protocol.stretchtag, '.mat'];
result_filename = ['Rush-Larsen/vw-results/vw-results_N', num2str(N),...
   '_t', tissue_string,...
   '_s', num2str(sample),...
   '_s2c',num2str(s2cellnum), protocol.drugtag, protocol.ascendtag, protocol.rjtag, protocol.stretchtag, '.mat'];

% check to see if ICs have been determined prior 
if isfile(result_filename)
   disp(['Already done tissue ', tissue_string,...
      ', sample ', num2str(sample),...
      ', s2c ', num2str(s2cellnum)])
   if minerva % Mount Sinai HPC, disregard
      blank = [];
      save(result_filename, "blank")
   end
   return
end

if ~isfile(ic_filename) % if no ICs gathered yet for this condition

   if isfile([ic_filename_partial, '_crash', protocol.drugtag, protocol.ascendtag, protocol.rjtag, protocol.stretchtag, '.mat'])
      disp(['Tissue ', tissue_string,...
         ', sample ', num2str(sample),...
         ' previously crashed'])
      return
   end
   %%
   disp('first pace at 1Hz for 21 beats')
   protocol.n_stim = 21;
   protocol.stim_delay = 10 ; % ms
   protocol.stim_dur = 1 ; % ms
   protocol.stim_amp = 210 ; %A/F
   protocol.stim_freq = 1; % Hz

   % AP is divided into two chunks, the first (350ms) with small dt when 
   % state variables are changing rapidly, and the second with larger dt 
   % when the state variables are changing more slowly. This compromise 
   % improves efficiency of the Crank-Nicolson method which cannot use a 
   % fully variable timestep
   protocol.slow_dur = 350; % ms
   protocol.dt_small = 5e-4; % set such that dV is less than 1 for all time steps
   protocol.dt_big = 2e-3; % ms
   protocol.stim_delay = 100 ; % ms
   
   % run simple cable simulation for 21 beats
   [t, ~, ~, state_all_last_beat] = f_tomek_cable_CNRL_CVTS(N, sample, tissue, protocol);
   
   % High heterogeneity cables occasionally crash due to numerical
   % instabilities. Below are some steps to address this issue. No cables
   % ultimately used in our study required this handling.
   if size(state_all_last_beat, 3) == 1
      % if crash once, try again with smaller time step
      disp(['Try tissue ', tissue_string, ', sample ', num2str(sample),' again (1)'])
      protocol.dt_small = 5e-4;
         protocol.dt_big = 5e-4;
         [t, ~, ~, state_all_last_beat] = f_tomek_cable_CNRL_CVTS(N, sample, tissue, protocol);

      if size(state_all_last_beat, 3) == 1
         % if crash again, try modifying the stimulus (longer duration,
         % lower amplitude, applied to CM #1-3.
         disp(['Try tissue ', tissue_string, ', sample ', num2str(sample),' again (2)'])
         protocol.stim_delay = 100 ; % ms
         protocol.stim_dur = 3 ; % ms
         protocol.stim_amp = 80 ;
         protocol.stim_cell = 1:3;
         [t, ~, ~, state_all_last_beat] = f_tomek_cable_CNRL_CVTS(N, sample, tissue, protocol);

         if size(state_all_last_beat, 3) == 1
            % if a third crash occurs, conclude the attempt
            warning(['could not pace tissue ', tissue_string,...
               ', sample ', num2str(sample),' to steady state'])
            t_s1=t;
            state_all_s1=state_all_last_beat;
            save([ic_filename_partial, '_crash', protocol.drugtag, protocol.ascendtag, '.mat'], 't_s1', 'state_all_s1')
            return
         end
      end
   end

   t_s1 = t - t(1);
   state_all_s1 = state_all_last_beat;
   save(ic_filename, 't_s1', 'state_all_s1') % save ICs

   clear state_all_last_beat state_all_s1 t t_s1 % clear workspace
end

%% Once ICs have been established, we can begin prpearing to test S2:
initstepmax = 1.00001 * (vw(4) - vw(3));
% step sizes for the bifurcation algorithm
steps = [163.84, 81.92, 40.96, 20.48, 10.24, 5.12, 2.56, 1.28, 0.64, 0.32, 0.16, 0.08, 0.04, 0.02, 0.01];%, 0.005, 0.0025, 0.00125];
% only test steps smaller than the initial window
steps = steps(steps <= initstepmax);

%% Save conductances relevant to the current simulation if not done already
load('Populations/minicombo.mat', 'sigmas');
sigma = sigmas(tissue);

if ~isfield(protocol, 'stretch')
   stretchtag = '';
   protocol.stretch=[];
else
   stretchtag = protocol.stretchtag;
end

if ~isfile(['Populations/sample', num2str(sample),'_vw_conds_c', num2str(N),'_s', num2str(sigma),'_t', tissue_string, stretchtag,'.mat'])
    % all conductances for all populations are saved in minicombo.mat. This
    % section extracts the relevant conductances for this simulation.
    cond_tensor = load('Populations/minicombo.mat', 'cond_tensor');
   if tissue > 1
      cond_tensor = cond_tensor.cond_tensor(:,:,tissue); 
   else
      if sample > 1
         cond_tensor = cond_tensor.cond_tensor(:,:,uniform_variant);
      else
         cond_tensor = cond_tensor.cond_tensor(:,:,tissue);
      end
   end
   conds_selection = load(['Populations/sample',num2str(sample),'_conds_selection.mat'], 'conds_selection');
   conds_selection = conds_selection.conds_selection;
   conds_all = cond_tensor(conds_selection(1:N), :)'; % dimensions: ncells x nconds
   cond_tensor = conds_all;
   if tissue == 1
      cond_tensor = repmat(cond_tensor(:,1), 1, N); % make uniformly the first cell
   end

   % insert a uniform stretch in the middle of the cable if required, shave
   % off ends to maintain cable length
   if ~isempty(protocol.stretch)
      repCM = repmat(cond_tensor(:,50),1,length(protocol.stretch));
      nToTrim = (length(protocol.stretch)-1)/2;
      upstreamHet = cond_tensor(:,nToTrim+1:49);
      downstreamHet = cond_tensor(:,51:end-nToTrim);
      cond_tensor = [upstreamHet,repCM,downstreamHet];
   end
   
   % save conductances relevant to the current simulation for easy access
   save(['Populations/sample',num2str(sample),'_vw_conds_c', num2str(N),'_s', num2str(sigma),'_t', tissue_string, stretchtag, '.mat'], 'cond_tensor', 'sigmas')
   clear cond_tensor conds_all conds_selection
end

%% prepare placeholders for stored results
pre = [];
post = [];
mid = [];
pre_result = [];
post_result = [];
mid_result = struct('mid_intervals', []);

for st = 1:length(steps) % iterate through steps
   ready = [0, 0];
   while ~(ready(1) && ready(2)) % only advance to smaller step size once pre- and post- timepoints are found
      % initialize intervals for this iteration of the algorithm
      vw_ints = vw(3):steps(st):vw(4)*1.0000001;
      % remove intervals that have been decided in prior loops
      if ~isempty(pre_result)
         vw_ints = vw_ints(2:end) ;
      end
      if ~isempty(post_result)
         vw_ints = vw_ints(1:end-1) ;
      end
      if ~isempty(mid)
         vw_ints = [vw_ints(vw_ints < 0.9999999999*mid(1)), vw_ints(vw_ints > 1.000000001*mid(end))];
      end
      vw_full = [vw(1:2), vw_ints];

      %% Run simulations
      protocol.n_stim = 1;
      protocol.stim_dur = 1 ; % ms
      protocol.stim_amp = s2a ; 
      protocol.stim_freq = 1; % Hz
      protocol.slow_dur = 1; % ms
      protocol.dt_small = 1e-3  ; % timestep constant for VW determination
      protocol.dt_big = 1e-3;

      nvw = length(vw_ints); % how many S2 timepoints are tested on the current loop
      for ints = 1:nvw
         protocol.s2_time = vw_ints(ints);
            
         % abandon attempt if cannot find pre-VW state by 50ms or post-VW 
         % state by 700ms after S1
         if protocol.s2_time > 700 || protocol.s2_time < 50
            vw_failure = 1 ;
            break 
         end
         
         % run simulation
         [t, V_all, ~, ~, ij_data] = f_tomek_cable_CNRL_CVTS(N, sample, tissue, protocol, vw_full);
         
         %% Classify output by VW state
         window = 0;
         
         % to identify which state (pre-/mid-/post-VW) the cable is in,
         % look at whether propagation has reached the ends of the cable
         % (first/last 20%)
         last20pct = floor(0.8*length(t));

         % the threshold for identifying propagation in a region is 50mV
         % above the minimum voltage seen in that region
         vw_threshold_upstream = min(min(V_all(:,1:s2cellnum))) + 50;
         vw_threshold_downstream = min(min(V_all(:,s2cellnum:end))) + 50;

         if max(max(V_all(last20pct:end,1:floor(N/4)))) < vw_threshold_upstream ... 
               && max(max(V_all(last20pct:end,floor(3*N/4):N))) < vw_threshold_downstream
            % propagation has reached neither upstream nor downstream end
            window = 1; % pre-vw
            pre = vw_ints(ints);

         elseif max(max(V_all(last20pct:end,1:floor(N/4)))) > vw_threshold_upstream ...
               && max(max(V_all(last20pct:end,floor(3*N/4):N))) > vw_threshold_downstream
             % propagation has reached both upstream and downstream ends
            window = 3; % post-vw
            post = vw_ints(ints);

         elseif max(max(V_all(last20pct:end,1:floor(N/4)))) > vw_threshold_upstream ...
               && max(max(V_all(last20pct:end,floor(3*N/4):N))) < vw_threshold_downstream
            % propagation has reached upstream end but not downstream end
             window = 2; % mid-vw
            mid = vw_ints(ints);

         elseif max(max(V_all(last20pct:end,1:floor(N/4)))) < vw_threshold_upstream ...
               && max(max(V_all(last20pct:end,floor(3*N/4):N))) > vw_threshold_downstream
            % propagation has reached downstream end but not upstream end
             window = -2; % mid-rvs-vw
            mid = vw_ints(ints);
         else
            warning('window status not identified')
         end

         % display progress (optional):
         disp([tissue, sample, s2cellnum, steps(st), protocol.s2_time, window])

         %% store result for pre-vw if there is one
         if window == 1
            % if our run has produced a novel pre-vw result:
            pre_result = struct(...
               't', t,...
               'V_all', V_all,...
               'new', 1, 'pre', pre,...
               'ij_data', ij_data);
            ready(1) = 1;
            pre = protocol.s2_time;
         else
            % if a pre-vw result exists, but not new in this iteration:
            if ~isempty(pre_result)
             ready(1) = 1;
               pre_result.new = 0;
            end
         end

         %% store result from post-vw if there is one
         if window == 3
            % if our run has produced a novel post-vw result:
            post_result = struct(...
               't', t,...
               'V_all', V_all,...
               'new', 1, 'post', post,...
               'ij_data', ij_data);
            ready(2) = 1;
            post = protocol.s2_time;
         else
            % if a post-vw result exists, but not new in this iteration:
            if ~isempty(post_result)
               ready(2) = 1;
               post_result.new = 0;
            end
         end

         %% store result from mid-vw if there is one
         if abs(window) == 2
            % if our run has produced a novel mid-vw result:
            if ~isempty(mid)

               % if there was not previously a mid result
               if isempty(mid_result.mid_intervals)
                  mid_result = struct(...
                     't_open', t,...
                     'V_open', V_all,...
                     't_close', t,...
                     'V_close', V_all,...
                     'new', 1, 'mid_intervals', mid, 'type','null',...
                     'ij_data_open', ij_data, 'ij_data_close', ij_data);

                  % or if there was previously one mid result
               elseif length(mid_result.mid_intervals) == 1
                  % window extends later
                  if mid > mid_result.mid_intervals
                     mid_result.t_close = t;
                     mid_result.V_close = V_all;
                     mid_result.new = 1;
                     mid_result.mid_intervals = [mid_result.mid_intervals, mid];
                     mid_result.ij_data_close = ij_data;

                     % or window extends earlier
                  else
                     mid_result.t_open = t;
                     mid_result.V_open = V_all;
                     mid_result.new = 1;
                     mid_result.mid_intervals = [mid, mid_result.mid_intervals];
                     mid_result.ij_data_open = ij_data;
                  end

                  % or if there were previously two or more mid results
               elseif length(mid_result.mid_intervals) == 2
                  % window extends earlier
                  if mid < mid_result.mid_intervals(1)
                     mid_result.t_open = t;
                     mid_result.V_open = V_all;
                     mid_result.new = 1;
                     mid_result.mid_intervals(1) = mid;
                     mid_result.ij_data_open = ij_data;

                     % or window extends later
                  else
                     mid_result.t_close = t;
                     mid_result.V_close = V_all;
                     mid_result.new = 1;
                     mid_result.mid_intervals(2) = mid;
                     mid_result.ij_data_close = ij_data;
                  end
               end

               % save the type of vulnerable window; "normal" upstream
               % propagating (VW) or "reverse" downstream propagating
               % (RVW)
               if window == 2
                  mid_result.type = 'VW';
               elseif window == -2
                  mid_result.type = 'RVW';
               end

               % or if this run has NOT produced a novel mid-vw result
            else
               if ~isempty(mid_result)
                  mid_result.new = 0;
               end
            end
            mid = mid_result.mid_intervals;
         end

      end
      if vw_failure
         break
      end
      %% set new bounds for next loop
      pre_next = pre;
      post_next = post;
      if isempty(pre)
         pre_next = vw_full(3) - steps(st);
      elseif isempty(post)
         post_next = vw_full(end) + steps(st);
      end

      vw(3:4) = [pre_next, post_next];
   end
   if vw_failure
      break
   end
end

%%
vw_duration = 0;
if ~isempty(mid)
   vw_duration = mid(end)-mid(1);
end
disp('final result:')
if vw_failure
   disp('failed to elicit s2 AP')
else
   disp([sample, rj, s2a, vw_duration, pre, post])
end

%% trim results to reduce size of results file
[~,I] = max(pre_result.V_all(:,s2cellnum));
cutoff = 40/protocol.dt_small + I;
pre_result.t = pre_result.t(1:cutoff);
pre_result.V_all = pre_result.V_all(1:cutoff,:);

if isfield(mid_result, 'V_open')
   [~,I1] = max(mid_result.V_open(:,1));
   [~,I2] = max(mid_result.V_open(:,end));
   I = max([I1, max(I2)]);
   cutoff = 10/protocol.dt_small + I;
   if cutoff < length(mid_result.t_open)
      mid_result.t_open = mid_result.t_open(1:cutoff);
      mid_result.V_open = mid_result.V_open(1:cutoff,:);
   end

   [~,I1] = max(mid_result.V_close(:,1));
   [~,I2] = max(mid_result.V_close(:,end));
   I = max([I1, max(I2)]);
   cutoff = 10/protocol.dt_small + I;
   if cutoff < length(mid_result.t_close)
      mid_result.t_close = mid_result.t_close(1:cutoff);
      mid_result.V_close = mid_result.V_close(1:cutoff,:);
   end
end

[~,I1] = max(post_result.V_all(:,1));
[~,I2] = max(post_result.V_all(:,end-10:end));
I = max([I1, max(I2)]);
cutoff = 10/protocol.dt_small + I;
if cutoff < length(post_result.t)
   post_result.t = post_result.t(1:cutoff);
   post_result.V_all = post_result.V_all(1:cutoff,:);
end

%% save results

if minerva % Mount Sinai HPC, disregard
   save(['/sc/arion/scratch/pullit02/vw-results/vw-results_N', num2str(N),...
      '_t', tissue_string,...
      '_s', num2str(sample),...
      '_s2c',num2str(s2cellnum), protocol.drugtag,  protocol.rjtag,'.mat'],'protocol', 'N', 's2cellnum', 'sample', 'tissue', 'sigma','vw',...
      'pre', 'pre_result','mid', 'mid_result', 'post', 'post_result', 'vw_duration', 'vw_failure')
   blank = [];
   save(result_filename, 'blank')
else
   save(result_filename, 'protocol', 'N', 's2cellnum', 'sample', 'tissue', 'sigma','vw',...
      'pre', 'pre_result','mid', 'mid_result', 'post', 'post_result', 'vw_duration')
end
end
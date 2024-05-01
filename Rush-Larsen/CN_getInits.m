function [v, nai, nass, ki, kss, cai, cass, cansr, cajsr, m, hp, h, j, jp, mL,...
   hL, hLp, a, iF, iS, ap, iFp, iSp, d, ff, fs, fcaf, fcas, jca, nca, nca_i,...
   ffp, fcafp, xs1, xs2, Jrel_np, CaMKt, ikr_c0, ikr_c1, ikr_c2, ikr_o, ikr_i,...
   Jrel_p] = CN_getInits(filename, N, partway, protocol)

%% if starting at all baseline:
if isequal(filename, 'Rush-Larsen/ICs/CN_inits.mat')

   load(filename, 'v', 'nai', 'nass', 'ki', 'kss', 'cai', 'cass', 'cansr',...
      'cajsr', 'm', 'hp', 'h', 'j', 'jp', 'mL', 'hL', 'hLp', 'a', 'iF', 'iS',...
      'ap', 'iFp', 'iSp', 'd', 'ff', 'fs', 'fcaf', 'fcas', 'jca', 'nca', 'nca_i',...
      'ffp', 'fcafp', 'xs1', 'xs2', 'Jrel_np', 'CaMKt', 'ikr_c0', 'ikr_c1',...
      'ikr_c2', 'ikr_o', 'ikr_i', 'Jrel_p')

   v     = v     * ones(N,1) ;
   nai   = nai   * ones(N,1) ;
   nass  = nass  * ones(N,1) ;
   ki    = ki    * ones(N,1) ;
   kss   = kss   * ones(N,1) ;
   cai   = cai   * ones(N,1) ;
   cass  = cass  * ones(N,1) ;
   cansr = cansr * ones(N,1) ;
   cajsr = cajsr * ones(N,1) ;
   m     = m     * ones(N,1) ;
   hp    = hp    * ones(N,1) ;
   h     = h     * ones(N,1) ;
   j     = j     * ones(N,1) ;
   jp    = jp    * ones(N,1) ;
   mL    = mL    * ones(N,1) ;
   hL    = hL    * ones(N,1) ;
   hLp   = hLp   * ones(N,1) ;
   a     = a     * ones(N,1) ;
   iF    = iF    * ones(N,1) ;
   iS    = iS    * ones(N,1) ;
   ap    = ap    * ones(N,1) ;
   iFp   = iFp   * ones(N,1) ;
   iSp   = iSp   * ones(N,1) ;
   d     = d     * ones(N,1) ;
   ff    = ff    * ones(N,1) ;
   fs    = fs    * ones(N,1) ;
   fcaf  = fcaf  * ones(N,1) ;
   fcas  = fcas  * ones(N,1) ;
   jca   = jca   * ones(N,1) ;
   nca   = nca   * ones(N,1) ;
   nca_i = nca_i * ones(N,1) ;
   ffp   = ffp   * ones(N,1) ;
   fcafp = fcafp * ones(N,1) ;
   xs1   = xs1   * ones(N,1) ;
   xs2   = xs2   * ones(N,1) ;
   Jrel_np = Jrel_np * ones(N,1) ;
   CaMKt = CaMKt * ones(N,1) ;
   ikr_c0 = ikr_c0    * ones(N,1) ;
   ikr_c1 = ikr_c1 * ones(N,1) ;
   ikr_c2 = ikr_c2  * ones(N,1) ;
   ikr_o = ikr_o    * ones(N,1) ;
   ikr_i = ikr_i   * ones(N,1) ;
   Jrel_p = Jrel_p * ones(N,1) ;

   %% if starting at pre_ICs:
elseif regexp(filename, regexptranslate('wildcard', 'Rush-Larsen/ICs/preics*.mat'))
   v     = ones(N,1) ;
   nai   = ones(N,1) ;
   nass  = ones(N,1) ;
   ki    = ones(N,1) ;
   kss   = ones(N,1) ;
   cai   = ones(N,1) ;
   cass  = ones(N,1) ;
   cansr = ones(N,1) ;
   cajsr = ones(N,1) ;
   m     = ones(N,1) ;
   hp    = ones(N,1) ;
   h     = ones(N,1) ;
   j     = ones(N,1) ;
   jp    = ones(N,1) ;
   mL    = ones(N,1) ;
   hL    = ones(N,1) ;
   hLp   = ones(N,1) ;
   a     = ones(N,1) ;
   iF    = ones(N,1) ;
   iS    = ones(N,1) ;
   ap    = ones(N,1) ;
   iFp   = ones(N,1) ;
   iSp   = ones(N,1) ;
   d     = ones(N,1) ;
   ff    = ones(N,1) ;
   fs    = ones(N,1) ;
   fcaf  = ones(N,1) ;
   fcas  = ones(N,1) ;
   jca   = ones(N,1) ;
   nca   = ones(N,1) ;
   nca_i = ones(N,1) ;
   ffp   = ones(N,1) ;
   fcafp = ones(N,1) ;
   xs1   = ones(N,1) ;
   xs2   = ones(N,1) ;
   Jrel_np = ones(N,1) ;
   CaMKt = ones(N,1) ;
   ikr_c0 = ones(N,1) ;
   ikr_c1 = ones(N,1) ;
   ikr_c2 = ones(N,1) ;
   ikr_o = ones(N,1) ;
   ikr_i = ones(N,1) ;
   Jrel_p = ones(N,1) ;

   %load(filename, 'tissue_string', 'sample', 'tissue');
   tissue_string = protocol.tissue_string;
   sample = protocol.sample;
   tissue = protocol.tissue;

   tissue_full = str2double(tissue_string);
   uniform_variant = floor(10.00000001*(tissue_full-tissue));

   %% set up & save conductances for this tissue sample
   load('Populations/minicombo.mat', 'sigmas');
   sigma = sigmas(tissue);


   if ~isfield(protocol, 'stretch')
      stretchtag = '';
      protocol.stretch=[];
   else
      stretchtag = protocol.stretchtag;
   end


   if ~isfile(['Populations/sample', num2str(sample),'_vw_conds_c', num2str(N),'_s', num2str(sigma),'_t', tissue_string, stretchtag, '.mat'])
      cond_tensor = load('Populations/minicombo.mat', 'cond_tensor');
      if tissue > 1
         cond_tensor = cond_tensor.cond_tensor(:,:,tissue); % saved as a struct
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

if ~isempty(protocol.stretch)
repCM=repmat(cond_tensor(:,50),1,length(protocol.stretch));
nToTrim=(length(protocol.stretch)-1)/2;
upstreamHet = [cond_tensor(:,nToTrim+1:49)];
downstreamHet = [cond_tensor(:,51:end-nToTrim)];
cond_tensor = [upstreamHet,repCM,downstreamHet];
end




      save(['Populations/sample',num2str(sample),'_vw_conds_c', num2str(N),'_s', num2str(sigma),'_t', tissue_string, stretchtag,'.mat'], 'cond_tensor', 'sigmas')
      clear conds_all conds_selection
   else
      load(['Populations/sample',num2str(sample),'_vw_conds_c', num2str(N),'_s', num2str(sigma),'_t', tissue_string, stretchtag,'.mat'], 'cond_tensor')
   end

conductances = cond_tensor';
conductances = conductances .* protocol.drugs;
if protocol.ascend > 0
    conductances = sortrows(conductances, protocol.ascend);
elseif protocol.ascend < 0
    conductances = sortrows(conductances, protocol.ascend*-1, 'descend');
end
cond_tensor = conductances';

   init_Y = f_getInit_Y('endo');
   CM = ToRORd19_fast(init_Y, 0, 0) ;
   CM.ODEModel = @f_ToRORd19 ;
   setUpPacingProtocol(CM,40,...
      80,...
      999,...
      1);
   if tissue > 1
      for i = 1:N
         scaleConductances(CM, cond_tensor(:,i)', CM.conductanceNames)

         odeSolver(CM);
         CM.conductances.scaling = CM.conductances.baseline;

         v(i)     = CM.state.Y(end,1) ;
         nai(i)   = CM.state.Y(end,2) ;
         nass(i)  = CM.state.Y(end,3) ;
         ki(i)    = CM.state.Y(end,4) ;
         kss(i)   = CM.state.Y(end,5) ;
         cai(i)   = CM.state.Y(end,6) ;
         cass(i)  = CM.state.Y(end,7) ;
         cansr(i) = CM.state.Y(end,8) ;
         cajsr(i) = CM.state.Y(end,9) ;
         m(i)     = CM.state.Y(end,10) ;
         hp(i)    = CM.state.Y(end,11) ;
         h(i)     = CM.state.Y(end,12) ;
         j(i)     = CM.state.Y(end,13) ;
         jp(i)    = CM.state.Y(end,14) ;
         mL(i)    = CM.state.Y(end,15) ;
         hL(i)    = CM.state.Y(end,16) ;
         hLp(i)   = CM.state.Y(end,17) ;
         a(i)     = CM.state.Y(end,18) ;
         iF(i)    = CM.state.Y(end,19) ;
         iS(i)    = CM.state.Y(end,20) ;
         ap(i)    = CM.state.Y(end,21) ;
         iFp(i)   = CM.state.Y(end,22) ;
         iSp(i)   = CM.state.Y(end,23) ;
         d(i)     = CM.state.Y(end,24) ;
         ff(i)    = CM.state.Y(end,25) ;
         fs(i)    = CM.state.Y(end,26) ;
         fcaf(i)  = CM.state.Y(end,27) ;
         fcas(i)  = CM.state.Y(end,28) ;
         jca(i)   = CM.state.Y(end,29) ;
         nca(i)   = CM.state.Y(end,30) ;
         nca_i(i) = CM.state.Y(end,31) ;
         ffp(i)   = CM.state.Y(end,32) ;
         fcafp(i) = CM.state.Y(end,33) ;
         xs1(i)   = CM.state.Y(end,34) ;
         xs2(i)   = CM.state.Y(end,35) ;
         Jrel_np(i)= CM.state.Y(end,36) ;
         CaMKt(i)  = CM.state.Y(end,37) ;
         ikr_c0(i) = CM.state.Y(end,38) ;
         ikr_c1(i) = CM.state.Y(end,39) ;
         ikr_c2(i) = CM.state.Y(end,40) ;
         ikr_o(i)  = CM.state.Y(end,41) ;
         ikr_i(i)  = CM.state.Y(end,42) ;
         Jrel_p(i) = CM.state.Y(end,43) ;
      end
   else
      scaleConductances(CM, cond_tensor(:,1)', CM.conductanceNames)

         odeSolver(CM);
         CM.conductances.scaling = CM.conductances.baseline;

         v     =  repmat(CM.state.Y(end,1),N,1) ;
         nai   =  repmat(CM.state.Y(end,2),N,1) ;
         nass  =  repmat(CM.state.Y(end,3),N,1) ;
         ki    =  repmat(CM.state.Y(end,4),N,1) ;
         kss   =  repmat(CM.state.Y(end,5),N,1) ;
         cai   =  repmat(CM.state.Y(end,6),N,1) ;
         cass =  repmat(CM.state.Y(end,7),N,1) ;
         cansr =  repmat(CM.state.Y(end,8),N,1) ;
         cajsr =  repmat(CM.state.Y(end,9),N,1) ;
         m     =  repmat(CM.state.Y(end,10),N,1) ;
         hp    =  repmat(CM.state.Y(end,11),N,1) ;
         h     =  repmat(CM.state.Y(end,12),N,1) ;
         j     =  repmat(CM.state.Y(end,13),N,1) ;
         jp    =  repmat(CM.state.Y(end,14),N,1) ;
         mL    =  repmat(CM.state.Y(end,15),N,1) ;
         hL    =  repmat(CM.state.Y(end,16),N,1) ;
         hLp   =  repmat(CM.state.Y(end,17),N,1) ;
         a     =  repmat(CM.state.Y(end,18),N,1) ;
         iF    =  repmat(CM.state.Y(end,19),N,1) ;
         iS    =  repmat(CM.state.Y(end,20),N,1) ;
         ap    =  repmat(CM.state.Y(end,21),N,1) ;
         iFp   =  repmat(CM.state.Y(end,22),N,1) ;
         iSp   =  repmat(CM.state.Y(end,23),N,1) ;
         d     =  repmat(CM.state.Y(end,24),N,1) ;
         ff    =  repmat(CM.state.Y(end,25),N,1) ;
         fs    =  repmat(CM.state.Y(end,26),N,1) ;
         fcaf  =  repmat(CM.state.Y(end,27),N,1) ;
         fcas  =  repmat(CM.state.Y(end,28),N,1) ;
         jca   =  repmat(CM.state.Y(end,29),N,1) ;
         nca   =  repmat(CM.state.Y(end,30),N,1) ;
         nca_i =  repmat(CM.state.Y(end,31),N,1) ;
         ffp   =  repmat(CM.state.Y(end,32),N,1) ;
         fcafp =  repmat(CM.state.Y(end,33),N,1) ;
         xs1   =  repmat(CM.state.Y(end,34),N,1) ;
         xs2   =  repmat(CM.state.Y(end,35),N,1) ;
         Jrel_np = repmat(CM.state.Y(end,36),N,1) ;
         CaMKt  =  repmat(CM.state.Y(end,37),N,1) ;
         ikr_c0 =  repmat(CM.state.Y(end,38),N,1) ;
         ikr_c1 =  repmat(CM.state.Y(end,39),N,1) ;
         ikr_c2 =  repmat(CM.state.Y(end,40),N,1) ;
         ikr_o  =  repmat(CM.state.Y(end,41),N,1) ;
         ikr_i  =  repmat(CM.state.Y(end,42),N,1) ;
         Jrel_p =  repmat(CM.state.Y(end,43),N,1) ;
   end

   if isfile(filename)
   delete(filename) %delete pre-ICs file
   end

elseif regexp(filename, regexptranslate('wildcard', 'Rush-Larsen/ICs/ics*beat21*.mat'))
   load(filename, 'state_all_s1', 't_s1')
   if ~isempty(partway)
      dex = find(abs(t_s1-partway) < 1e-5);
   else
      dex = length(t_s1);
   end

   v  = squeeze(state_all_s1(:,1,dex));
   nai = squeeze(state_all_s1(:,2,dex));
   nass = squeeze(state_all_s1(:,3,dex));
   ki = squeeze(state_all_s1(:,4,dex));
   kss = squeeze(state_all_s1(:,5,dex));
   cai = squeeze(state_all_s1(:,6,dex));
   cass = squeeze(state_all_s1(:,7,dex));
   cansr = squeeze(state_all_s1(:,8,dex));
   cajsr = squeeze(state_all_s1(:,9,dex));
   m = squeeze(state_all_s1(:,10,dex));
   hp = squeeze(state_all_s1(:,11,dex));
   h = squeeze(state_all_s1(:,12,dex));
   j = squeeze(state_all_s1(:,13,dex));
   jp = squeeze(state_all_s1(:,14,dex));
   mL = squeeze(state_all_s1(:,15,dex));
   hL = squeeze(state_all_s1(:,16,dex));
   hLp = squeeze(state_all_s1(:,17,dex));
   a = squeeze(state_all_s1(:,18,dex));
   iF = squeeze(state_all_s1(:,19,dex));
   iS = squeeze(state_all_s1(:,20,dex));
   ap = squeeze(state_all_s1(:,21,dex));
   iFp = squeeze(state_all_s1(:,22,dex));
   iSp = squeeze(state_all_s1(:,23,dex));
   d = squeeze(state_all_s1(:,24,dex));
   ff = squeeze(state_all_s1(:,25,dex));
   fs = squeeze(state_all_s1(:,26,dex));
   fcaf = squeeze(state_all_s1(:,27,dex));
   fcas = squeeze(state_all_s1(:,28,dex));
   jca = squeeze(state_all_s1(:,29,dex));
   nca = squeeze(state_all_s1(:,30,dex));
   nca_i = squeeze(state_all_s1(:,31,dex));
   ffp = squeeze(state_all_s1(:,32,dex));
   fcafp = squeeze(state_all_s1(:,33,dex));
   xs1 = squeeze(state_all_s1(:,34,dex));
   xs2 = squeeze(state_all_s1(:,35,dex));
   Jrel_np = squeeze(state_all_s1(:,36,dex));
   CaMKt = squeeze(state_all_s1(:,37,dex));
   ikr_c0 = squeeze(state_all_s1(:,38,dex));
   ikr_c1 = squeeze(state_all_s1(:,39,dex));
   ikr_c2 = squeeze(state_all_s1(:,40,dex));
   ikr_o = squeeze(state_all_s1(:,41,dex));
   ikr_i = squeeze(state_all_s1(:,42,dex));
   Jrel_p = squeeze(state_all_s1(:,43,dex));

   if length(v) ~= N
      disp([filename, " LENGTH V = ", num2str(length(v)),"; N = ",num2str(N)])
      save([filename(1:end-4),'IC_BUG.mat'])
      %error('saved initial conditions have wrong cable length')
   end
end
end

classdef ToRORd19_fast < AbsCM
    properties
        % defined in @f_ToRORd19_loadBaselineProperties
        universals; celltype; name ; ODEModel ; isCellAlive ; cellID ; geometry ; state ; YNames ; YUnits ;
        currentNames ; currentUnits ; parameterNames ; parameterUnits ;
        parameters ; conductanceNames ; conductances ; protocol ;
        x_conductances ; y_allevents ; x_allevents ; event_number ; isPaced ;
    end

    methods
        function tord19cm = ToRORd19_fast(init_Y, celltype, cellID)
            [universals, name, model, isCellAlive, geometry, state, YNames, YUnits,...
                currentNames, currentUnits, parameterNames, parameterUnits,...
                parameters, conductanceNames, conductances, protocol, isPaced]...
                = f_ToRORd19_loadBaselineProperties(init_Y) ;
            tord19cm.universals = universals ;
            tord19cm.celltype = celltype ;
            tord19cm.ODEModel = model ;
            tord19cm.isCellAlive = isCellAlive ;
            tord19cm.cellID = cellID ;
            tord19cm.geometry = geometry ;
            tord19cm.state = state ;
            tord19cm.parameterNames = parameterNames ;
            tord19cm.parameters = parameters ;
            tord19cm.conductanceNames = conductanceNames ;
            tord19cm.conductances = conductances ;
            tord19cm.protocol = protocol ;
            tord19cm.isPaced = isPaced ;
        end

        function saveX_conductance(tord19cm, x_conductance)
            tord19cm.x_conductances = x_conductance ;
        end

        function scaleParameters(tord19cm, scaleFactors, names)
            [~, indA, indB] = intersect(names, tord19cm.parameterNames, 'stable') ;
            tord19cm.parameters.scaling(indB) = tord19cm.parameters.scaling(indB) .* scaleFactors(indA);
        end

        function scaleConductances(tord19cm, scaleFactors, names)
            [~, indA, indB] = intersect(names, tord19cm.conductanceNames, 'stable');
            tord19cm.conductances.scaling(indB) = tord19cm.conductances.scaling(indB) .* scaleFactors(indA);
        end

        function setEnvironment(tord19cm, T, nao, cao, ko)
            tord19cm.universals.T = T ;
            tord19cm.universals.nao = nao ;
            tord19cm.universals.cao = cao ;
            tord19cm.universals.ko = ko ;
        end


        function setUpPacingProtocol(tord19cm, amplitudes, numPulses, precedingTime, pulseDurations)
            [tord19cm.protocol.intervalTimes, tord19cm.protocol.stimulus]...
                = f_setUpPacingProtocol(amplitudes, numPulses, precedingTime, pulseDurations) ;

            tord19cm.protocol.amplitudes = amplitudes ;
            tord19cm.protocol.numPulses = numPulses ;
            tord19cm.protocol.precedingTime = precedingTime ;
            tord19cm.protocol.pulseDuration = pulseDurations ;

            tord19cm.protocol.phaseLengths = numPulses.*(precedingTime + pulseDurations) ; % for reference
            tord19cm.protocol.phaseLengths(1) = 0.2 + (numPulses(1) - 1) * precedingTime(1) + numPulses(1) * pulseDurations(1) ; % first pulse doesn't have full preceding time
            tord19cm.protocol.phaseLengths(end) = tord19cm.protocol.phaseLengths(end) + 1 ; % add extra second at the end to capture last AP
            tord19cm.protocol.frequencies = 1000./(precedingTime + pulseDurations) ; % for reference

            tord19cm.protocol.totalTime = tord19cm.protocol.intervalTimes(end,2) ;
            tord19cm.protocol.nPhases = length(amplitudes) ;
            tord19cm.protocol.phaseLengths = numPulses.*(precedingTime + pulseDurations)/1000 ;
            tord19cm.protocol.intervalTimes = tord19cm.protocol.intervalTimes * 1000; % change to ms
        end

        function odeSolver(tord19cm)

            if isempty(tord19cm.protocol.intervalTimes)
                error('Error: Please set up pacing protocol first')
            else
                options = odeset('MaxStep', 1, 'InitialStep', 2e-2);% solver settings

                % Allocate space initially to save time. Extra zeros will be
                % removed later to save space.
                initialrows = 500000 ; % enough for 3 beats

                tord19cm.state.Y = zeros(initialrows, length(tord19cm.state.init_Y)) ;
                tord19cm.state.t = zeros(initialrows, 1) ; % first row at T0, no need to update

                tord19cm.state.Y(1,:) = tord19cm.state.init_Y ; % first row, input values (generally steady state)
                input_values = tord19cm.state.init_Y ; % initialize, changes in for-loop
                index = 1 ; % initialize, changes in for-loop

                stimIndices = find(tord19cm.protocol.stimulus > 0);
                i_stimToKeep = 0;

                % convert properties stored as structs to cells to pass to
                % ode15s
                universals_arr = cell2mat(struct2cell(tord19cm.universals));
                geometry_arr = cell2mat(struct2cell(tord19cm.geometry));
                parameters_arr = tord19cm.parameters.baseline.*tord19cm.parameters.scaling ;
                conductances_arr = tord19cm.conductances.baseline.*tord19cm.conductances.scaling ;

                timer = 0;
                for i = 1:length(tord19cm.protocol.stimulus)
                    tic;
                    if tord19cm.isCellAlive == 1
                        input_values(end) = tord19cm.protocol.stimulus(i) ; % last input value is the stimulus amplitude
                        [t_current, Y_current] = ode15s(...
                            tord19cm.ODEModel,...
                            tord19cm.protocol.intervalTimes(i,:),...
                            input_values,...
                            options,...
                            tord19cm.celltype,...
                            universals_arr,...
                            geometry_arr,...
                            parameters_arr,...
                            conductances_arr);

                        current_length = length(t_current) ; % length of this interval in the stimulus protocol

                    end

                    % cells may have been killed off above so we check again for
                    % survival
                    if tord19cm.isCellAlive == 1
                        if i >= i_stimToKeep
                            tord19cm.state.t(index:index+current_length-1) = t_current(1:end) ; % save outcome of this run
                            tord19cm.state.Y(index:index+current_length-1,:) = Y_current(1:end,:) ; % save outcome of this run
                            index = index + current_length-1; % update the index only if storing all or if last 3 beats
                        end
                        input_values = tord19cm.state.Y(index, :) ; % set input values for next interval to end of current interval
                    end
                    timer = timer + toc ;

                end

                if ~isempty(tord19cm.state.t) && tord19cm.isCellAlive == 1
                    % remove extra zeros to save space
                    zeroInds = find(~tord19cm.state.t) ; % returns indices of all zeroes
                    firstExtraZero = zeroInds(2) ; % first zero is T0
                    tord19cm.state.t = tord19cm.state.t(1:firstExtraZero-1) ;
                    tord19cm.state.Y = tord19cm.state.Y(1:firstExtraZero-1,:) ;

                    tord19cm.state.t = tord19cm.state.t(2:end);
                    tord19cm.state.Y = tord19cm.state.Y(2:end, :);
                    tord19cm.protocol.last3stimTimes = tord19cm.protocol.intervalTimes(stimIndices(end-2:end),1) - tord19cm.state.t(1);
                    tord19cm.state.t = tord19cm.state.t - tord19cm.state.t(1); % set t0 to 0
                end
            end
        end

        function getCurrents(tord19cm)
            if isempty(tord19cm.state.t)
                error('Error: Please run the odeSolver on this ToRORd19 CM first')
            else
                tord19cm.state.currents = zeros(length(tord19cm.state.t), length(tord19cm.currentNames));
                numLoops = length(tord19cm.state.t) ;

                % convert properties stored as structs to cells to pass to
                % ode15s
                universals_arr = cell2mat(struct2cell(tord19cm.universals));
                geometry_arr = cell2mat(struct2cell(tord19cm.geometry));
                parameters_arr = tord19cm.parameters.baseline.*tord19cm.parameters.scaling ;
                conductances_arr = tord19cm.conductances.baseline.*tord19cm.conductances.scaling ;

                for i= 1:numLoops
                    [~, ij_data]    = f_ToRORd19(tord19cm.state.t(i), tord19cm.state.Y(i,:),...
                        tord19cm.celltype,...
                        universals_arr,...
                        geometry_arr,...
                        parameters_arr,...
                        conductances_arr); % ~ says only return data, not dY output
                    tord19cm.state.currents(i,1:end-1) = ij_data;
                end
                tord19cm.state.currents(:,end) = tord19cm.state.currents(:,8) + tord19cm.state.currents(:,9) ; % INaCa = INaCa_i + INaCa_ss
            end
        end
    end
end
function [instCV2080, meanCV2080, wideCV2080] = calculateCV(t, V)%t:its x 1; V:its x ncells
%%
[~, ncells] = size(V); % get number of cells total
pctl20 = floor(ncells*0.2)-1; % which cell is at 20% -1
pctl80 = ceil(ncells*0.8)+1; % which cell is at 80% +1
posV = zeros(1,pctl80-pctl20+1);
for i = pctl20:pctl80
    posVs = find(V(:,i) > 0); % find index where upstroke reaches cell i (V>0)
    if isempty(posVs)
    posV(i-pctl20+1) = 1;
    else
    posV(i-pctl20+1) = posVs(1);
    end
end
%%
t_posV1 = t(posV); % time where upstroke reaches each cell
t_posV2 = [t_posV1(3:end);0;0]; % time where upstroke reaches 2 cells later
t_posVdiff = t_posV2(1:end-2)-t_posV1(1:end-2); % inst time for upstroke to move 2 cells

instCV2080 = 0.02./(t_posVdiff./1000); % instant CV for each cell distance (cm/s)
meanCV2080 = mean(instCV2080); % mean CV of all instant CVs
wideCV2080 = 0.01*(pctl80-pctl20)./((t_posV1(end)-t_posV1(1))./1000); % overall CV based on time from cell 20%-1 to cell 80%+1

end

    
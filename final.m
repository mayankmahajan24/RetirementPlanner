% - - - - - PART 1 (Scenario Generation - - - - - - %
run1 = 1; %Generate returns
run2 = 0; %CVaR frontier
run3 = 0; %Bogle
run4 = 0; %Trigger
run5 = 0;
run6 = 0;
run7 = 0;
run8 = 0;
%weights = [0 0.02 0.98];
%[cap, varisk, cvar] = simulate(all_returns, num_years_71, weights, invest_pct, [], [], []);
%hist(cap);
%mean(cap)
%sum(cap >= 7.5*10^6)

if run1
    'Running Part 1'
    all_returns = load('all_returns.csv');
    weights = [0.0 0.4 0.6];

    ann_normal = [ .007 .035 .122];
    ann_crash= [.007 .115 -.18];

    Qnormal = [ .0001/246	0	0;
    0	4.73936E-06	-5.26929E-06;
    0	-5.26929E-06	9.6077E-05] * 246;

    Qcrash = [.0001/246	0 0;
    0	1.10541E-05	-2.63882E-05;
    0	-2.63882E-05	0.000470094] * 246;

    Fnormal = chol(Qnormal);
    Fcrash = chol(Qcrash);

    TM = [1 0; .5 .5];
    T = 50;

    num_scenarios = 1000; % Number of scenarios
    cum_capital = zeros(num_scenarios,50);
    num_years = 50;
    n_assets = 3;
    %all_returns = zeros(num_scenarios * num_years, n_assets); 

    for x = 1:num_scenarios
        CDF = cumsum(TM,2); % cumulative distribution function.
        i = 1; % index of the initial shock.
        current_sal = 125000;
        current_invest = current_sal * 0.16; %initial capital
        prev_capital = 0;

        for t = 1:T % simulation for T periods.     
            if i == 1 %normal
                ret = ann_normal + randn(1 ,n_assets) * Fnormal; 
            elseif i == 2 %crash
                ret = ann_crash + randn(1 ,n_assets) * Fcrash;
            else print 'what!'
            end

            ret = max(ret, -1);

            %all_returns((t-1)*num_scenarios + x,:) = ret;
            
            annual = all_returns((t-1)*num_scenarios + x,:)*transpose(weights);
            current_invest = current_sal * 0.16;
            cum_capital(x,t) = (prev_capital + current_invest) * (1 + annual);
            prev_capital = cum_capital(x,t);
            current_sal = current_sal * 1.04; %salary increase every year
            
            ug(t) = rand(1); % create pseudo-random numbers from the uniform dist.
            j = find(CDF(i,:)>=ug(t),1,'First');
            
            if j ~= 1 & j ~= 2
                print 'WHAT'
            end
            i = j;
        end      
    end





    r67 = cum_capital(:,36);
    r71 = cum_capital(:,40);

    inds = [];
    for i = 1:1000
       if ~isempty(find(cum_capital(i,:) > 7.5 * 10^6, 1))
           inds = [inds find(cum_capital(i,:) > 7.5 * 10^6, 1)];
       end
    end

    var67 = -quantile(r67, 0.05)
    on = r67 < quantile(r67, 0.05);
    ES67 = - sum(on .* r67)/sum(on)

    var71 = -quantile(r71, 0.05)
    on = r71 < quantile(r71, 0.05);
    ES71 = - sum(on .* r71)/sum(on)

    %MAX DRAWDOWN
    [MaxDD, MaxDDIndex] = maxdrawdown(transpose(cum_capital));
    %maximumDrawDown = max(MaxDD)
    %theIndex = find(MaxDD == max(MaxDD), 1)
    %st = max(MaxDDIndex(1,theIndex)-2,1)
    %en = min(MaxDDIndex(2,theIndex)+2,50)
    %plot(st:en, cum_capital(theIndex, st:en));
    %cum_capital(theIndex, st:en)
    sorted = sort(MaxDD);
    avMDD = mean(sorted(951:1000))

end

% ------------------PART 2 (CVar Frontier) ----------------- %

if run2
    point_space = zeros(1000000,3);
    in = 1;
    for cash = 0:0.01:1
        for bonds = 0:0.01:1
            for stock = 0:0.01:1
                if cash + bonds + stock <= 1.0
                    point_space(in,:) = [cash bonds stock];
                    in = in + 1;
                end
            end
        end
    end

    point_space = point_space(1:in-1,:);
    n = 1000;
    datanum = 3;
    r_exp = 10^6 : 0.5*10^6 : 8 * 10^6;
    %r_exp = 10^6;
    num_scenarios = 1000;
    n_assets = 3;
    invest_pct = 0.16;

    opt_port = zeros(length(r_exp),datanum); % Save the optimal portfolios
    vol_port = zeros(length(r_exp),2);       % Save the corresponding volatility
    cvar_port = zeros(length(r_exp),2);      % Save the corresponding VaR and cVaR
    min_mean = zeros(length(r_exp),1);
    num_years = 36;

    opt_port_71 = zeros(length(r_exp),datanum); % Save the optimal portfolios
    vol_port_71 = zeros(length(r_exp),2);       % Save the corresponding volatility
    cvar_port_71 = zeros(length(r_exp),2);      % Save the corresponding VaR and cVaR
    min_mean_71 = zeros(length(r_exp),1);
    num_years_71 = 40;


    %return scenarios are in all_returns

    for j = 1:length(r_exp)                  % Optimize w.r.t. each point on efficient frontier
        % cvx_begin quiet
        %     variables weights(datanum) var1 slack(n) % thres is the VaR (alpha in equation(25).)
        %     expressions cap(num_scenarios) r_all(num_scenarios*num_years, n_assets) current_sal prev_capital(num_scenarios)
        %     minimize var1 + sum(slack)/(n*.02)
        %     subject to
        %     %Set up some variables
        %     current_invest = 0; 
        %     current_sal = 125000;
        %     for i = 1:num_years %for each year.
        %         current_invest = current_sal * 0.16; %Update new savings
        %         returns = all_returns(num_scenarios*(i-1) + 1 : num_scenarios*i, :);
        %         annual = returns * weights; %Get overall growth
                
        %         cap = (cap + current_invest) .* (1 + annual); %For each scenario, evolve the capital
        %         current_sal = current_sal * 1.04; %salary increase every year
        %     end
        %     sum(weights) <= 1                     % Budget constraint
        %     weights >= 0
        %     sum(cap) >= num_scenarios * r_exp(j)     % Expected return constraint
        %     slack >= 0                               % slack's are auxiliary variables;
        %     cap + var1 + slack >= 0      % they stand for max(...,0)
        % cvx_end
        j

        min_cvar = 999999999;
        min_cvar_71 = 999999999;
        for pt = 1:length(point_space)
            weights = point_space(pt,:);
            %Set up some variables  
            
            [cap, valatrisk, cvar] = simulate(all_returns, num_years, weights, invest_pct, [], [], []);

            if (mean(cap) >= r_exp(j) & cvar <= min_cvar)            
                opt_port(j,:) = weights';
                vol_port(j,2) = sqrt(var(cap,1));
                vol_port(j,1) = vol_port(j,2);     % Annualized volatility
                cvar_port(j,1) = cvar;                  % CVaR
                cvar_port(j,2) = valatrisk;  % VaR, can also use "thres"
                min_cvar = cvar;
                min_mean(j) = mean(cap);
            end

            [cap, valatrisk_71, cvar_71] = simulate(all_returns, num_years_71 - num_years, weights, invest_pct, cap, num_years+1, []);

            
            if (mean(cap) >= r_exp(j) & cvar_71 <= min_cvar_71)            
                opt_port_71(j,:) = weights';
                vol_port_71(j,2) = sqrt(var(cap,1));
                vol_port_71(j,1) = vol_port_71(j,2);     % Annualized volatility
                cvar_port_71(j,1) = cvar_71;                  % CVaR
                cvar_port_71(j,2) = valatrisk_71;  % VaR, can also use "thres"
                min_cvar_71 = cvar_71;
                min_mean_71(j) = mean(cap);
            end
        end

        cvar_port(j,1)
        min_mean(j)

    end

    %plot(-cvar_port(6:11,1), min_mean(6:11));

    %plot(-cvar_port_71(6:11,1), min_mean_71(6:11));
end

%----------------PART 3 BOGLE RULE--------------------------%

if run3
    stock_pct = 0:0.1:1;
    means67 = zeros(11,1);
    cvars67 = zeros(11,1);

    means71 = zeros(11,1);
    cvars71 = zeros(11,1);

    for i = 1:length(stock_pct)
        weights = [0 1-stock_pct(i) stock_pct(i)];

        %evaluate mean and cvar
        [cap, varisk67, cvars67(i)] = simulate(all_returns, num_years, weights, invest_pct, [], [], []);
        means67(i) = mean(cap);
        [cap, varisk71, cvars71(i)] = simulate(all_returns, num_years_71 - num_years, weights, invest_pct, cap, num_years+1, []);
        means71(i) = mean(cap);
    end

    %  plot(cvars67, means67);
    %  for i = 1:length(stock_pct)
    %      text(cvars67(i), means67(i), [num2str(stock_pct(i)*100) '/' num2str(100-stock_pct(i)*100)]);
    %  end
    % % plot(cvars71, means71);
    % for i = 1:length(stock_pct)
    %     text(cvars71(i), means71(i), [num2str(stock_pct(i)*100) '/' num2str(100-stock_pct(i)*100)]);
    % end

    %Test bogles rule 1000 scenarios

    weights = [0 0.31 0.69];

    [cap, varisk67bogle, cvar67bogle] = simulate(all_returns, num_years, weights, invest_pct, [], [], 'Bogle');
    mean67bogle = mean(cap);

    weights = [0 0.31+(num_years)*0.01 0.69-(num_years)*0.01];

    [cap, varisk71bogle, cvar71bogle] = simulate(all_returns, num_years_71 - num_years, weights, invest_pct, cap, num_years+1, 'Bogle');
    mean71bogle = mean(cap);

end

%----------------PART 4 TRIGGERS --------------------------%



if run4
    FMTmeans = zeros(length(opt_port_71), 1);
    FMTvars = zeros(length(opt_port_71), 1);
    FMTcvars = zeros(length(opt_port_71), 1);
    FMTgoals = zeros(length(opt_port_71), 3);

    for i = 1:length(opt_port_71)
        weights = opt_port_71(i,:);
        [cap_all, FMTvars(i), FMTcvars(i), FMTgoals(i,:)] = simulateTrigger(all_returns, num_years_71, weights, invest_pct, [], [], []);
        FMTmeans(i) = mean(cap_all(:,num_years_71));
    end

    weights = [0 0.31 0.69];
    [BT_cap_all, BTvar, BTcvar, BTgoals] = simulateTrigger(all_returns, num_years_71, weights, invest_pct, [], [], 'Bogle');
    BTmeans = mean(BT_cap_all(:,num_years_71));

end

% ----------------PART 5 SAVINGS RATE ----------------------%

if run5
    invest_range = 0.16:0.01:0.30;
    s_means = zeros(length(invest_range), 1);
    s_vars = zeros(length(invest_range), 1);
    s_cvars = zeros(length(invest_range), 1);
    s_goals = zeros(length(invest_range), 3);
    weights = [0 0.38 0.62];
    for i = 1:length(invest_range)
        %weights is highest risk fix mix
        [cap_all, s_vars(i), s_cvars(i), s_goals(i,:)] = simulateTrigger(all_returns, num_years_71, weights, invest_range(i), [], [], []);
        s_means(i) = mean(cap_all(:, num_years_71));
    end

    plot(invest_range, s_goals(:,1));
end


% -----------------PART 6 MAYANKS RULE -------------------------%


if run6
    nvest_range = 0.16:0.01:0.30;
    c_means = zeros(length(invest_range), 1);
    c_vars = zeros(length(invest_range), 1);
    c_cvars = zeros(length(invest_range), 1);
    c_goals = zeros(length(invest_range), 3);
    weights = [0 0 1];
    for i = 1:length(invest_range)
        %weights is highest risk fix mix
        [cap_all, c_vars(i), c_cvars(i), c_goals(i,:)] = simulateTrigger(all_returns, num_years_71, weights, invest_range(i), [], [], 'Custom');
        c_means(i) = mean(cap_all(:, num_years_71));
    end
    plot(invest_range, c_goals(:,1));
end

if run7
    invest_range = 0.16:0.01:0.30;
    ce_means = zeros(length(invest_range), 1);
    ce_vars = zeros(length(invest_range), 1);
    ce_cvars = zeros(length(invest_range), 1);
    ce_goals = zeros(length(invest_range), 3);
    weights = [0 0 1];
    for i = 1:length(invest_range)
        %weights is highest risk fix mix
        [cap_all, ce_vars(i), ce_cvars(i), ce_goals(i,:)] = simulateTrigger(all_returns, num_years_71, weights, invest_range(i), [], [], 'CustomExtra');
        ce_means(i) = mean(cap_all(:, num_years_71));
    end
    plot(invest_range, ce_goals(:,1));
end

if run8
    [cap_all, p_varisk, p_cvar, p_goals] = simulateTrigger(all_returns, num_years_71+4, weights, 0.16, [], [], 'Custom');
    [cap_all, pe_varisk, pe_cvar, pe_goals] = simulateTrigger(all_returns, num_years_71+4, weights, 0.16, [], [], 'CustomExtra');
end




function [cap, varisk, cvar] = simulate(all_returns, num_years, weights, invest_pct, cap_init, year_init, simType)
    current_sal = 125000;
    num_scenarios = 1000;
    cap = zeros(num_scenarios, 1);
    cap_all = zeros(num_scenarios, 50);

    year_start = 1;
    goals = zeros(num_scenarios, 1);
    goalsCount = 1;

    if nargin >= 5 & ~isempty(cap_init)
        cap = cap_init;
    end

    if nargin >= 6 & ~isempty(year_init)
        year_start = year_init;
    end

    for i = year_start:(num_years+year_start-1) %for each year.
        current_invest = current_sal * invest_pct; %Update new savings
        returns = all_returns(num_scenarios*(i-1) + 1 : num_scenarios*i, :);
        annual = returns * transpose(weights); %Get overall growth
        cap = (cap + current_invest) .* (1 + annual); %For each scenario, evolve the capital
        current_sal = current_sal * 1.04; %salary increase every year

        if ~isempty(simType) & string(simType) == string('Bogle')
            weights(2) = weights(2) + 0.01;
            weights(3) = weights(3) - 0.01;
        end
        cap_all(:,i) = cap;
    end

    varisk = -quantile(cap,0.05);
    cvar = -1* sum((cap <= -varisk) .* cap) / sum(cap <= -varisk);
end


function [cap_all, varisk, cvar, FMTgoals] = simulateTrigger(all_returns, num_years, weights, invest_pct, cap_init, year_init, simType)
    current_sal = 125000;
    num_scenarios = 1000;
    cap_all = zeros(num_scenarios, 50);
    cap = zeros(num_scenarios, 1);
    year_start = 1;
    end_goal = 7.5 * 10^6;

    if nargin >= 5 & ~isempty(cap_init)
        cap = cap_init;
    end

    if nargin >= 6 & ~isempty(year_init)
        year_start = year_init;
    end

    for i = year_start:(num_years+year_start-1) %for each year.

        if ~isempty(simType) & string(simType) == string('CustomExtra') & 31 + i == 55
            invest_pct = invest_pct + 0.06;
        end

        current_invest = current_sal * invest_pct; %Update new savings
        returns = all_returns(num_scenarios*(i-1) + 1 : num_scenarios*i, :);
        
        if ~isempty(simType) & (string(simType) == string('Custom') | string(simType) == string('CustomExtra'))
            annual = zeros(num_scenarios, 1); 
            for x = 1:num_scenarios
                weights = [1 cap(x)/end_goal (1-cap(x)/end_goal)];
                annual(x) = returns(x,:) * transpose(weights);
            end
        else
            annual = returns * transpose(weights); %Get overall growth
        end
        %Trigger stuff
        reached_goal = (cap >= 7.5 *10^6);
        annual(reached_goal) = returns(reached_goal, 1); %zeros(sum(goal),1);

        cap = (cap + current_invest) .* (1 + annual); %For each scenario, evolve the capital
        current_sal = current_sal * 1.04; %salary increase every year

        if ~isempty(simType) & string(simType) == string('Bogle')
            weights(2) = weights(2) + 0.01;
            weights(3) = weights(3) - 0.01;
        end


        cap_all(:,i) = cap;
    end

    varisk = -quantile(cap,0.05);
    cvar = -1* sum((cap <= -varisk) .* cap) / sum(cap <= -varisk);

    triggerInds = [];

    for i = 1:1000
       if ~isempty(find(cap_all(i,:) >= 7.5 * 10^6, 1))
           triggerInds = [triggerInds find(cap_all(i,:) >= 7.5 * 10^6, 1)];
       end
    end

    FMTgoals = [length(triggerInds)/length(cap) mean(triggerInds) std(triggerInds)];

end



    


if ~exist('clear_flag', 'var'), clear_flag = 1; end

if usejava('desktop') && clear_flag
   clear; 
else
    ver;
end
close all;

%--------------------------------------------------------------------------
% simulation setup
%--------------------------------------------------------------------------

% file with model
respath='./';
outpath='./Results/';
if ~exist('resfile_list','var')
    resfile_list={'res_20191112_agg','res_20191112_om25'};
        
end

for f=1:length(resfile_list)
    
    resfile=resfile_list{f};
    
    load([respath,resfile,'.mat']);
    
    % Initial Economy Config
    varlist={'simseries','statevec','indexmap','varnames'};
    load(['sim_',resfile],varlist{:});
    
    % set starting point
    start_ini=5;
    start_shock=[0,3,4]; %[0,3,4];
    statevec=statevec(2:end);
    startvals=mean(simseries(statevec==start_ini,:));
    N_vars=length(startvals);
    % states or variables computed at market prices?
    
    
    % number of periods and burn-in
    N_shock=length(start_shock);
    N_runs=10000;
    NT_sim=25;
    NT_ini=0;
    NT_sim=NT_sim+1;
    
    % cluster endogenous states
    if mobj.Params.consolidatedLender
        envarind=4:7;
    else
        envarind=4:8;
    end
    maxclust=10;
    enstatemat=simseries(statevec==start_ini,envarind);
    cindex=clusterdata(enstatemat,'criterion','distance','maxclust',maxclust,'linkage','weighted');
    sttot=size(enstatemat,1);
    startptmat=[];
    for c=1:maxclust
        cfrac=sum(cindex==c)/sttot;
        thismean=mean(enstatemat(cindex==c,:),1);
        disp([num2str(c),': ',num2str(thismean)]);
        thisc=repmat(thismean,floor(N_runs*cfrac),1);
        startptmat=[startptmat; thisc];
    end
    if size(startptmat,1)<N_runs
        thismean=mean(enstatemat);
        startptmat=[startptmat; repmat(thismean,N_runs-size(startptmat,1),1)];
    end
    
    
    % report levels or grwoth rates for output variables
    reportLevels=1;
    
    % compute Euler equation error?
    compEEErr=0;
    % make graphs grayscale
    grayscale=0;
    
    % Compute term premium (slow!)
    term_premium=0;
    
    % output table file
    outfile=['GR_',resfile];
    
    varnames_store = varnames;
    
    simseries_median = cell(N_shock,1);
    simseries_mean = cell(N_shock,1);
    simseries_std = cell(N_shock,1);
    
    simseries_diff_median = cell(N_shock,1);
    simseries_diff_mean = cell(N_shock,1);
    simseries_diff_std = cell(N_shock,1);
    clear simseries;

    open_parpool;
    
    
    for s=1:N_shock
        
        disp(['Shock ',num2str(s),' of ',num2str(N_shock)]);
        tens_simseries = zeros(NT_sim,N_vars,N_runs);
        
        % compute entry of random number matrix that sets first state
        % deterministically to start_shock
        if start_shock(s)>0
            transprob=cumsum(mobj.Exogenv.mtrans(start_ini,:));
            shock_prob=transprob(start_shock(s));
            if start_shock(s)>1
                shock_prob_minus=transprob(start_shock(s)-1);
            else
                shock_prob_minus=0;
            end
            rvar_next=(shock_prob+shock_prob_minus)/2;
        end
        
        % Create shock matrix
        rng(1);
        %shmatfull = rand(NT_sim*N_runs,1);
        shmatfull = lhsdesign(N_runs,NT_sim);
        
        SDFmat=zeros(NT_sim,mobj.Exogenv.exnpt);
        
        fprintf([repmat('.',1,100) '\n\n']);
        
        parfor n=1:N_runs
            %--------------------------------------------------------------------------
            % start simulation
            %--------------------------------------------------------------------------
            %fprintf('Run %d - Start \n',n);
            % simulate
            shmat = shmatfull(n,:)';
            if start_shock(s)>0
                shmat(1)=rvar_next;
            end
            
            startpt=struct;
            startpt.AB=startptmat(n,1);
            startpt.MB=startptmat(n,2);
            startpt.KREO=startptmat(n,3);
            if ~mobj.Params.consolidatedLender
                startpt.WI=startptmat(n,4);
                startpt.BG=startptmat(n,5);
            else
                startpt.BG=startptmat(n,4);                
            end
            startpt=orderfields(startpt,mobj.En_names);
            startpt_vec=model.DSGEModel.structToVec(startpt)';
            startpt_vec=[start_ini,startpt_vec];
            
            [simseries,varnames,~,~,~]=mobj.simulate(NT_sim,NT_ini,startpt_vec,compEEErr,shmat);
            simseries_orig=simseries;
            varnames_orig=varnames;
            statevec = simseries(:,1);
            
            [simseries, varnames] = mobj.computeSimulationMoments(simseries,varnames);
            
            nvars = length(varnames);
            tens_simseries(:,:,n) = [startvals; simseries];
            if mod(n,N_runs/100)==0
                %disp([num2str(n),'/',num2str(N_runs),': ',num2str(round(1000*n/N_runs)/10),'% complete']);
                fprintf('\b|\n');
            end
            
        end
        fprintf('\n');
        varnames = varnames_store;
        nvars = length(varnames);
        
        % make HashMap with mapping of names to indices
        indexmap=java.util.HashMap;
        for i=1:nvars
            indexmap.put(varnames{i},i);
        end
        %     varst=zeros(length(startpt_vec)-1,1);
        %     for i=1:length(startpt_vec)-1
        %         varst(i)=indexmap.get(mobj.En_names{i});
        %     end
        
        %save(outfile,'tens_simseries','indexmap');
        
        simseries_median{s} = median(tens_simseries,3);
        simseries_mean{s} = mean(tens_simseries,3);
        simseries_std{s} = std(tens_simseries,[],3);
        
        if start_shock(s) > 0
            % If actual shock, difference and save
            tens_simseries_diff = tens_simseries - tens_simseries_0;
            simseries_diff_median{s} = median(tens_simseries_diff,3);
            simseries_diff_mean{s} = mean(tens_simseries_diff,3);
            simseries_diff_std{s} = std(tens_simseries_diff,[],3);
        else
            % If no shock, store for later differencing
            tens_simseries_0 = tens_simseries;
        end
        
    end
    
    save(outfile,'simseries_mean','simseries_median','simseries_std', ...
        'simseries_diff_mean', 'simseries_diff_median', 'simseries_diff_std', ...
        'indexmap','NT_sim','N_shock');
    
end

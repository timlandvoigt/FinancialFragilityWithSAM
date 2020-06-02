if ~exist('clear_flag', 'var'), clear_flag = 1; end

if usejava('desktop') && clear_flag
   clear; 
end
close all;

%--------------------------------------------------------------------------
% simulation setup
%--------------------------------------------------------------------------

% file with model
respath='./';
if ~exist('resfile','var')
    disp('resfile not defined. Opening default file instead.');
    resfile='res_20191112_om25'; 
end
load([respath,resfile,'.mat']);


% number of periods and burn-in
NT_sim=10000;
NT_ini=100;

% compute Euler equation error?
compEEErr=1;

% Winsorize
winsorize=0;
winsfile='';
cutoff=99.9;


% Force the creation of a sim_res file
force_output = 1;

% output table file
if ~exist('output_dir', 'var')
    output_dir = './Results/';
end

if ~exist(output_dir, 'dir')
    mkdir(output_dir)
end

%outpath='./Results_c/';
outstats_exog=[output_dir,'statsexog_',resfile,'.xls'];
errstats=[output_dir,'errstats_',resfile,'.xls'];
expdata=0;
outdata=[output_dir,'series_',resfile,'.csv'];
       
%--------------------------------------------------------------------------
% start simulation
%--------------------------------------------------------------------------

% set starting point
start_ex=3;
startpt=struct;
startpt.AB=stv{1}.State.AB;
startpt.MB=stv{1}.State.MB;
startpt.KREO=stv{1}.State.KREO;
if ~mobj.Params.consolidatedLender
    startpt.WI=stv{1}.State.WI;
end
startpt.BG=stv{1}.State.BG;
startpt=orderfields(startpt,mobj.En_names);
startpt_vec=model.DSGEModel.structToVec(startpt)';
startpt_vec=[start_ex,startpt_vec];

% simulate
[simseries,varnames,errmat,Wshtrans,SDFmat]=mobj.simulate(NT_sim+1,NT_ini,startpt_vec,compEEErr);
simseries_orig=simseries;
varnames_orig=varnames;
statevec = simseries(:,1);


[simseries, varnames] = mobj.computeSimulationMoments(simseries,varnames);
nvars = length(varnames);

% Create table object for easier access
simtable=array2table(simseries);
[~,ia,~]=unique(varnames);
simtable=simtable(:,ia);
simtable.Properties.VariableNames=varnames(ia);
dispnames=varnames(ia);

% make HashMap with mapping of names to indices
indexmap=java.util.HashMap;
for i=1:nvars
    indexmap.put(varnames{i},i);
end

% Check transition function errors
if compEEErr
    idx = sub2ind([NT_sim+1,mobj.Exogenv.exnpt],(1:NT_sim+1)',[statevec(2:end);1]);
    idx=idx(1:end-1);
    ABtrans=Wshtrans(:,1:mobj.Exogenv.exnpt);
    AB_err=simseries(:,indexmap.get('AB')) - ABtrans(idx); % column index: 6
    errmat = [errmat, [AB_err;0]];
    MBtrans=Wshtrans(:,mobj.Exogenv.exnpt+1:end);
    MB_err=simseries(:,indexmap.get('MB')) - MBtrans(idx);
    errmat = [errmat, [MB_err;0]];
    KREOtrans=Wshtrans(:,2*mobj.Exogenv.exnpt+1:end);
    KREO_err=simseries(:,indexmap.get('KREO')) -KREOtrans(idx);
    errmat = [errmat, [KREO_err;0]];
    WItrans=Wshtrans(:,3*mobj.Exogenv.exnpt+1:end);
    WI_err=simseries(:,indexmap.get('WI')) - WItrans(idx);
    errmat = [errmat, [WI_err;0]];
    BGtrans=Wshtrans(:,4*mobj.Exogenv.exnpt+1:end);
    BG_err=simseries(:,indexmap.get('BG')) - BGtrans(idx);
    errmat = [errmat, [BG_err;0]];
end

if winsorize
    if isempty(winsfile)
        prct=prctile(abs(errmat(2:end,end)),cutoff);
%        bad_idx = any(abs(errmat(2:end,end)) > repmat(prct,NT_sim,1) , 2);
        bad_idx = any(simseries(:,6) < repmat(.045,NT_sim,1) , 2);
    else
       load([respath,winsfile,'.mat'],'bad_idx');
    end
    simseries = simseries(~bad_idx,:);
    errmat = errmat([true;~bad_idx],:);
    NT_sim = size(simseries,1) + 1;
    save([respath,resfile,'.mat'],'bad_idx','-append');
end

%--------------------------------------------------------------------------
% calculate stats
%--------------------------------------------------------------------------
varst=zeros(length(startpt_vec)-1,1);
for i=1:length(startpt_vec)-1
    varst(i)=indexmap.get(mobj.En_names{i});
end
            
% state variable means in stationary distribution
stvstat=mean(simseries(:,varst));

% calculate business cycle stats
% first for all periods, then separately for low and high sigma_omega states
% 1. condition on exogenous states (smpsel_exog)
% 2. condition on endogenous states: GDP growth, default rate (smpsel_endog)
statsout_exog=cell(4,1);
statsout_endog=cell(4,1);

% first for all periods, then separately for low and high sigma_omega states
smpsel_exog={true(NT_sim-1,1), simseries(:,2)==mobj.Params.sig2_om(1) & simseries(:,1) >= mobj.Params.mu_Y, ...
                          simseries(:,2)==mobj.Params.sig2_om(1)  & simseries(:,1) < mobj.Params.mu_Y, ...
                          simseries(:,2)==mobj.Params.sig2_om(2)  & simseries(:,1) < mobj.Params.mu_Y, ...
                          simseries(:,2)==mobj.Params.sig2_om(1)};


gdp_idx = indexmap.get('Y');                     
gdpgr_idx = indexmap.get('Y_gr'); 
for j=1:numel(smpsel_exog)
    % subsample
    simtmp=simseries(smpsel_exog{j},:);
    statstmp=zeros(nvars,11);
    statstmp(:,1)=nanmean(simtmp)';
    statstmp(:,2)=nanstd(simtmp)';
    % contemp and first-order autocorrelations
    autocorrm=corrcoef([simtmp(2:end,:),simtmp(1:end-1,:)]);
    conm=autocorrm(1:nvars,1:nvars);
    lagm=autocorrm(nvars+1:end,1:nvars);
    % corr with shocks
    statstmp(:,3:4)=[conm(:,1),lagm(1,:)'];
    statstmp(:,5:6)=[conm(:,2),lagm(2,:)'];
    % corr with Y_gr
    statstmp(:,7:8)=[conm(:,gdp_idx),lagm(gdp_idx,:)'];
    statstmp(:,9:10)=[conm(:,gdpgr_idx),lagm(gdpgr_idx,:)'];
    % vector with fo autocorr
    statstmp(:,11)=diag(lagm);
    statsout_exog{j}=statstmp;
    
end


%--------------------------------------------------------------------------
% output
%--------------------------------------------------------------------------

% overview output for eyeball check against analytic st.st. values
% make one big structure with steady-state values
stvbig=model.HelperCollection.combineStructs({stv{1}.Sol,stv{1}.State,stv{1}.Add,stv{1}.statsout});

% output table
% make index vector
[displist,dispnames]=model.HelperCollection.makeListFromNames(indexmap,dispnames);
ndvars=length(displist);

disp(' ');
disp('Simulation steady state');

% overview output 
fprintf('Frequency (exog subsamples): ');
for j=1:numel(smpsel_exog)
    % select vars
    tabout_exog{j}=statsout_exog{j}(displist,:);
    fprintf('%f\t',sum(smpsel_exog{j}));
end
fprintf('\n');
disp('-------------');

for s=1:ndvars
    if isfield(stvbig,dispnames{s})
        ststval=stvbig.(dispnames{s});
    else
        ststval=0;
    end
    if numel(dispnames{s}) > 7
        fprintf('%d\t%4s\t\t\t\t%f |',displist(s),dispnames{s},ststval);
    else
        fprintf('%d\t%4s\t\t\t\t\t%f |',displist(s),dispnames{s},ststval);
    end
%     disp('Exog subsamples')
    for j=1:numel(smpsel_exog)
        fprintf('\t%f, %f |',tabout_exog{j}(s,1),tabout_exog{j}(s,2));
    end
    fprintf('\n');    
end

if compEEErr
    avg_err=mean(abs(errmat))';
    med_err=median(abs(errmat))';
    p75_err=prctile(abs(errmat),75)';
    p95_err=prctile(abs(errmat),95)';
    p99_err=prctile(abs(errmat),99)';
    p995_err=prctile(abs(errmat),99.5)';
    max_err=max(abs(errmat))';
    errtab=table(avg_err,med_err,p75_err,p95_err,p99_err,p995_err,max_err);
    errarr=table2array(errtab);
    disp(' ');
    disp('-----------------------------------------------');
    disp('Average and maximum Euler equation error');
    fprintf('Equ.no.\t\tAvg.\t\tMed.\t\tp75\t\t\tp95\t\t\tp99\t\t\tp99.5\t\tMax.\n');
    for s=1:length(avg_err)
        fprintf('%d\t\t\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n',s,errarr(s,1),errarr(s,2),errarr(s,3), ...
            errarr(s,4),errarr(s,5),errarr(s,6),errarr(s,7));
    end
    
    % plot EE error for these equations
    plotEE_pol=[3,4];
    plotEE_state=[0,0];
    for i=1:length(plotEE_pol)
        points=simseries(:,[6,7]);
        errvals=abs(errmat(1:end-1,plotEE_pol(i)));
        if plotEE_state(i)>0
            itmp=(statvec==plotEE_state(i));
            points=points(itmp,:);
            errvals=errvals(itmp,:);
        end
        model.HelperCollection.scatterPoints2D(points,errvals);
    end
    
end

% check grid bounds
state_range=4:8;
min_vec=min(simseries(:,state_range));
max_vec=max(simseries(:,state_range));
disp('State bounds:');
disp(mobj.Pfct.SSGrid.StateBounds(:,2:end));
disp('Simulation mins:');
disp(min_vec);
disp('Simulation max:');
disp(max_vec);


% write to file
values=struct2cell(mobj.Params);
paramout=cell2table(values,'RowNames',fieldnames(mobj.Params));
colnames={'mean','std','corrG','corrG_1','corrOm','corrOm_1','corrY','corrY_1','corrY_gr','corrY_gr_1','AC'};
for j=1:4
    tableout_exog=array2table(tabout_exog{j},'RowNames',dispnames,'VariableNames',colnames);
    writetable(tableout_exog,outstats_exog,'WriteRowNames',1,'FileType','spreadsheet','Sheet',j);
end
writetable(paramout,outstats_exog,'WriteRowNames',1,'FileType','spreadsheet','Sheet','params');
writetable(errtab,errstats,'FileType','spreadsheet');
if force_output
    params=mobj.Params;
    disp(['Saving simulation data to .mat file: ',['sim_',resfile,'.mat']]);
    save(['sim_',resfile,'.mat'],'simseries','displist','dispnames','errmat','tabout_exog','outstats_exog', ...
                'errstats','errtab','indexmap','NT_ini','NT_sim','smpsel_exog','statevec','statsout_exog','varnames','params');
end    

if expdata
    disp(' ');
    disp('Exporting simseries...');
    model.HelperCollection.tableExport(outdata,varnames,simseries);
end

% save model file with stationary state values
save([respath,resfile,'.mat'],'stvstat','-append');


if ~exist('clear_flag', 'var'), clear_flag = 1; end

if usejava('desktop') && clear_flag
   clear; 
end

close all;
respath='./';
outpath='./Results/';


% Benchmark vs. Agg Index

% if ~exist('econ1', 'var'), econ1='bench'; end
% if ~exist('econ2', 'var'), econ2='agg'; end
% if ~exist('colors', 'var'), colors = {'k-o', 'r-o'}; end
% if ~exist('label', 'var'), label = {'No Index', 'Agg. Index'}; end

% Benchmark vs. Local Only

econ1='bench';
econ2='om25';
colors = {'k-o', 'b-o'};
label = {'No Index', 'Local Index'};

% % Gov't Debt
% if ~exist('econ1', 'var'), econ1='bench'; end
% if ~exist('econ2', 'var'), econ2='aggfulltax'; end
% if ~exist('econ3', 'var'), econ3='aggtauL040'; end
% if ~exist('colors', 'var'), colors = {'k-o', 'r-o', 'b-o'}; end
% if ~exist('label', 'var'), label = {'No Index', 'Agg. Index', 'Agg + Tax Smoothing'}; end

% Benchmark vs. Agg Index

% if ~exist('econ1', 'var'), econ1='agg_i150'; end
% if ~exist('econ2', 'var'), econ2='aggtauL080_i150c100'; end
% if ~exist('colors', 'var'), colors = {'k-o', 'b-o'}; end
% if ~exist('label', 'var'), label = {'full tax', '\tau=0.8'}; end


if exist('econ3', 'var')
    N_economy = 3;
else
    N_economy=2; % numbers economies to plot on same graph
end


% econ2='xiom25indexonly';
%econ3='xi92benchgrid';
if ~exist('resfile','var')
    resfile=['res_20191112_finrec_',econ1,econ2];
    if ~exist('resfile1', 'var'), resfile1=['res_20191112_',econ1]; end
    if ~exist('resfile2', 'var'), resfile2=['res_20191112_',econ2]; end
    if N_economy == 3
        resfile = [resfile, econ3];
        if ~exist('resfile3', 'var'), resfile3 = ['res_20191112_',econ3]; end
    end
end
outfile=['GR_',resfile];
grayscale=0;
reportLevels=0;
batchMode=0;

load([respath,resfile1,'.mat']);
load([respath,'sim_',resfile1,'.mat']);
load([respath,'GR_',resfile1,'.mat']);

simseries_mean_econ{1} = simseries_mean{3}; % financial recession 
clear simseries_mean;

load([respath,resfile2,'.mat']);
load([respath,'sim_',resfile2,'.mat']);
load([respath,'GR_',resfile2,'.mat']);

simseries_mean_econ{2} = simseries_mean{3};
clear simseries_mean;

if N_economy == 3
    load([respath,resfile3,'.mat']);
    load([respath,'sim_',resfile3,'.mat']);
    load([respath,'GR_',resfile3,'.mat']);

    simseries_mean_econ{3} = simseries_mean{3};
    clear simseries_mean;
end

% load([respath,resfile3,'.mat']);
% load([respath,'sim_',resfile3,'.mat']);
% load([respath,'GR_',resfile3,'.mat']);
% 
% simseries_mean_econ{3} = simseries_mean{3};
% clear simseries_mean;

close all;

%% IRFS

outpath2=[outpath,outfile,'_'];
if ~exist('for_slides', 'var'), for_slides = 0; end
if for_slides
    outpath2 = [outpath2, 'slides_']; 
end

% file names for graphs (set to empty for no printing)
printfiles={[outpath2,'IRF1'],[outpath2,'IRF2'],[outpath2,'IRF3'],[outpath2,'IRF4'],[outpath2,'IRF5'],[outpath2,'IRF6']};        

brsel1=[indexmap.get('cB'),indexmap.get('cI'),indexmap.get('cS'),...
        indexmap.get('p'),indexmap.get('Drate_ZN'),indexmap.get('Drate_I')];
levelind1=ones(1,6);    
brsel2=[indexmap.get('Y'),indexmap.get('cB'),indexmap.get('cI'),...
         indexmap.get('cS'),indexmap.get('MB_g'),indexmap.get('BInext')];
levelind2=ones(1,6);
brsel3=[indexmap.get('Drate_ZN'),indexmap.get('Lspr_star'),indexmap.get('rD_real'),...
        indexmap.get('p'),indexmap.get('WI'),indexmap.get('Drate_I')];
levelind3=ones(1,6);
% brsel4=[indexmap.get('cB'),indexmap.get('Mstar'),indexmap.get('taxB'),indexmap.get('MB')];
% levelind4=ones(1,4);
brsel4=[indexmap.get('cB'), indexmap.get('netisspayB'),indexmap.get('issB'), ...
    indexmap.get('taxB'), indexmap.get('net_housing_B'),indexmap.get('net_rent_B')];
levelind4=ones(1,length(brsel4));
brsel5=[indexmap.get('cB'),indexmap.get('cI'),indexmap.get('rD_real'),...
        indexmap.get('p'),indexmap.get('Drate_ZN'),indexmap.get('Drate_I')];
levelind5=ones(1,6);   
brsel6=[indexmap.get('MB_g'), indexmap.get('fracMS'),indexmap.get('BInext'),indexmap.get('p'), ...
    indexmap.get('WI'), indexmap.get('cy'),indexmap.get('rD_real'),indexmap.get('pREO')];
levelind6=ones(1,8);


titles1={'Consumption B','Consumption I','Consumption D',...
         'House price','Loan defaults','Bank failures'}; 
titles2={'Output','Consumption B','Consumption I','Consumption S','Mortgage debt','Deposits'}; 
titles3={'Def. rate','Mortgage spread','Real riskfree',...
         'House price','Bank equity','Bank failures'};
% titles4={'Consumption B','M^*','Tax B'}; 
titles4 = {'Consumption B', 'Iss - Pay B', 'Issuance B', ...
         'Tax B', 'Net Housing B', 'Net Rent B'};
titles5={'Consumption B','Consumption I','Real riskfree',...
         'House price','Loan defaults','Bank failures'}; 
    
titles6 = {'Mortgage debt', 'MBS Savers (share)', 'Deposits','House price', ...
         'Bank equity', 'Conv. yield', 'Risk free real rate','REO house price'};


tvec=0:NT_sim-1;    
    
brsel_all=[brsel1,brsel2,brsel3,brsel4,brsel5,brsel6];
levelind_all=[levelind1, levelind2,levelind3, levelind4,levelind5,levelind6];
nvar=length(brsel_all);   
brseries_gr=zeros(N_economy, NT_sim, nvar);
for s=1:N_economy
    for v=1:nvar
        thisv=brsel_all(v);
        if levelind_all(v) 
            brseries_gr(s,:,v) = simseries_mean_econ{s}(:,thisv);
        else
            brseries_gr(s,:,v) = 100*(simseries_mean_econ{s}(:,thisv) ./ repmat(simseries_mean_econ{s}(1,thisv),NT_sim,1)-1);
        end
    end
end

% colors={'k-o','r-o','b-o'};
% label = {'no index.','aggreg. & 25% idios. index.'};
%label = {'no index.','aggreg. index.'};

% delta for mortgage debt
%brseries_gr(:,:,20)=brseries_gr(:,:,20)-brseries_gr(:,:,22);

if usejava('desktop')
    makeImpRes(brseries_gr(:,:,1:6),tvec,titles1,colors,[2,3],label,printfiles{1}, for_slides);   
    %makeImpRes(brseries_gr(:,:,7:12),tvec,titles2,colors,[2,3],label,printfiles{2}, for_slides);    
    %makeImpRes(brseries_gr(:,:,13:18),tvec,titles3,colors,[2,3],label,printfiles{3}, for_slides);    
    %makeImpRes(brseries_gr(:,:,19:24),tvec,titles4,colors,[2,3],label,printfiles{4}, for_slides);    
    %makeImpRes(brseries_gr(:,:,25:30),tvec,titles5,colors,[2,3],label,printfiles{5}, for_slides);    
    %makeImpRes(brseries_gr(:,:,31:end),tvec,titles6,colors,[2,4],label,printfiles{6}, for_slides);    
end

if batchMode
    close all;
    clear;
end


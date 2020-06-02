if ~exist('clear_flag', 'var'), clear_flag = 1; end

if usejava('desktop') && clear_flag
   clear; 
end

close all;
respath='./';
outpath='./Results/';
if ~exist('resfile','var')
    resfile='res_20191112_bench';
end
outfile=['GR_',resfile];
grayscale=0;
reportLevels=1;
batchMode=0;
relative_irf=1;

load([respath,resfile,'.mat']);
load([respath,'sim_',resfile,'.mat']);
load([respath,'GR_',resfile,'.mat']);

close all;

%% Make other Graphs
outpath=[outpath,outfile,'_'];
if relative_irf, outpath = [outpath, 'diff_']; end
tvec=0:NT_sim-1;

% file names for graphs (set to empty for no printing)
printfiles={[outpath,'IRF1'],[outpath,'IRF2'],[outpath,'IRF3'],[outpath,'IRF4']};        
% which variables
brsel1=[indexmap.get('Y'),indexmap.get('cB'),indexmap.get('cI'),...
         indexmap.get('cS'),indexmap.get('MB_g'),indexmap.get('BInext')];
brsel2=[indexmap.get('Drate_ZA'),indexmap.get('Lspr_star'),indexmap.get('rD_real'),...
        indexmap.get('p'),indexmap.get('WI'),indexmap.get('Drate_I')];
brsel3=[indexmap.get('cB'), indexmap.get('netisspayB'),indexmap.get('issB'), ...
    indexmap.get('taxB'), indexmap.get('net_housing_B'),indexmap.get('net_rent_B')];
brsel4=[indexmap.get('MB_g'), indexmap.get('MS_g'),indexmap.get('BInext'),indexmap.get('p'), ...
    indexmap.get('WI'), indexmap.get('cy'),indexmap.get('rD_real'),indexmap.get('pREO')];
% brsel3=[indexmap.get('LTV'),indexmap.get('Ifinlvg'),...
%         indexmap.get('Drate_I'),indexmap.get('bind_lamI')];

% How many shocks to plot (one less for relative)
if relative_irf
    N_shock_plot = N_shock - 1;
else
    N_shock_plot = N_shock;
end

brsel_all=[brsel1,brsel2,brsel3,brsel4];
nvar=length(brsel_all);   
brseries_gr=zeros(N_shock_plot, NT_sim, nvar);

for s=1:N_shock_plot
    if relative_irf
        if reportLevels==1
            brseries_gr(s,:,:) = simseries_diff_mean{s+1}(:,brsel_all);
        else
            brseries_gr(s,:,:) = 100*(simseries_diff_mean{s+1}(:,brsel_all) ./ repmat(simseries_diff_mean{s+1}(1,brsel_all),NT_sim,1)-1);
        end
    else
        if reportLevels==1
            brseries_gr(s,:,:) = simseries_mean{s}(:,brsel_all);
        else
            brseries_gr(s,:,:) = 100*(simseries_mean{s}(:,brsel_all) ./ repmat(simseries_mean{s}(1,brsel_all),NT_sim,1)-1);
        end
    end
end

colors={'b-o','r-o'};
if N_shock_plot==3
   colors = ['k-o',colors]; 
end
titles1={'Output','Consumption B','Consumption I','Consumption D','Mortgage debt','Deposits'}; 
titles2={'Def. rate','Mortgage spread','Risk free real rate',...
         'House price','Bank equity','Bank failures'};
titles3 = {'Consumption B', 'Iss - Pay B', 'Issuance B', ...
         'Tax B', 'Net Housing B', 'Net Rent B'};
titles4 = {'Mortgage debt', 'MBS Savers', 'Deposits','House price', ...
         'Bank equity', 'Conv. yield', 'Risk free real rate','REO house price'};
% titles3={'LTV','Interm. lvg', ...
%          'Interm. failures', 'bind_lamI'};  
     
if usejava('desktop')
    makeImpRes(brseries_gr(:,:,1:6),tvec,titles1,colors,[2,3],[],printfiles{1});   
    makeImpRes(brseries_gr(:,:,7:12),tvec,titles2,colors,[2,3],[],printfiles{2});
    makeImpRes(brseries_gr(:,:,13:18),tvec,titles3,colors,[2,3],[],printfiles{3});    
    makeImpRes(brseries_gr(:,:,19:end),tvec,titles4,colors,[2,4],[],printfiles{4});    
%    makeImpRes(brseries_gr(:,:,9:12),tvec,titles3,colors,[2,2],[],printfiles{3});  
%    makeImpRes(brseries_gr(:,:,13:16),tvec,titles4,colors,[2,2],[],printfiles{4});    
end


% if batchMode
%     close all;
%     clear;
% end

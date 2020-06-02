params=struct;
% mean-reverting AR(1) log(TFP)
params.sig_y=.008;
params.rho_y=.977;
params.mu_y=-.5*params.sig_y^2/(1-params.rho_y^2);
params.mu_Y=1;

% deterministic trend growth
params.g=0.0; 
params.mu_G=exp(params.g); 

% Omega transition matrix
Omprob_recess=[0.975, 0.025;
                0.075, 0.925];
params.Omprob_recess=Omprob_recess(:);    
Omprob_boom=[0.975, 0.025;
                0.075, 0.925];
params.Omprob_boom=Omprob_boom(:);      


% Pre-calibrated parameters
params.psi = 1.0000000000000000;
params.xihi = 0.21;
params.xilo = 0.15;
params.pi = 1.0057;
params.delta = 0.99565;
params.tau = 0.147;
params.tauL = 1.0;
params.phiK = 0.85;
params.nuK = 0.00616;
% note: steady state price from benchmark model
params.nuK_p = params.nuK * 8.8;
params.nuREO = 0.022;
params.nuREO_p = params.nuREO * 8.8;
params.SREO = 0.1666666666666667;
params.mu_kappa = 0.37;
params.sig_kappa = 0.152;
params.log_K_bar_total = 0;

% preferences
% borrowers
params.chiB = 0.343;
params.shareB = 0.47;
params.shareHB = 0.47;
params.sigmaB=2;
params.betaB = 0.95;
% intermediaries
params.chiI=.05;
params.shareI=.062;
params.shareHI=.062;
params.sigmaI=2;
params.betaI=params.betaB;
% savers
params.sigmaS=2;
params.betaS=.998;
params.xiL=0;

% housing and mortgage technology
params.maint_by_price=true;
params.sigom_lo=.2;
params.sigom_hi=.25;
params.rho_om=0.977;
params.rebate_maint=true;
params.psiS0=0.2;
params.psiS1=5;


% liquidity default
params.theta = 0.0;
params.etaB = 0.0;
params.liq_thresh = 1.0;
params.internalize_liq = 1.0;


% intermediation technology
params.phiI=0.93;
params.sigma_eps=0.065;
params.zeta=0.05;
%params.zeta=0.08;
params.noBR=0;

% indexation
params.iota_om=0;
params.iota_p=0;
params.pgr_max=Inf;
params.pgr_min=-Inf;   
params.indexAB=1;
params.indexMB=1;
params.indexTail=0;

% Quadrature 
[nodes, weights] = GaussHermite(11);
params.nodes=nodes;
params.weights=weights;

% model without separate intermediary sector
params.consolidatedLender=false;

% set KBbar 
params.KBbar = exp(params.log_K_bar_total) * params.shareHB;



% Define grids
startgrid.label = 'ini0';
startgrid.ABpts=linspace(0.033, 0.043, 6);
startgrid.MBpts=linspace(2.2, 2.85, 7);
startgrid.KREOpts=linspace(0, 0.05, 6);
startgrid.WIpts=[.01,.03,.05,.06,.08,.11,.13,.15,.19,.21,.23];
startgrid.BGpts=0;

benchgrid=struct('startgrid', startgrid);

agggrid0=startgrid;
agggrid0.ABpts=linspace(0.03, 0.045, 6);
agggrid0.MBpts=linspace(2.0, 2.65, 7);
agggrid0.WIpts=[-.06,-.03,-.01,.01,.03,.05,.06,.08,.11,.13,.15,.19,.20,.3,.35];
agggrid=struct('startgrid', agggrid0);

sh3grid0=startgrid;
sh3grid0.ABpts=linspace(0.028, 0.044, 6);
sh3grid0.MBpts=linspace(2.05, 2.8, 7);
sh3grid0.WIpts=[-.1,-.03,-.01,.01,.03,.05,.06,.08,.11,.13,.15,.19,.20,.3,.35];
sh3grid=struct('startgrid', sh3grid0);

reqgrid=startgrid;
reqgrid.ABpts=linspace(0.03, 0.045, 6);
reqgrid.MBpts=linspace(1.9, 2.5, 7);
reqgrid.WIpts=[-.03,-.01,.01,.03,.05,.06,.08,.11,.13,.15,.19,.20,.3,.35];
reqgrid.KREOpts=linspace(0.00, 0.03, 6);
reqgrid=struct('startgrid', reqgrid);



reggrid=startgrid;
reggrid.ABpts=linspace(0.03, 0.045, 6);
reggrid.MBpts=linspace(2.2, 2.85, 7);
reggrid.WIpts=[-.06,-.03,-.01,.01,.03,.05,.06,.08,.11,.13,.15,.19,.20,.3,.35];
reggrid.KREOpts=linspace(0.00, 0.03, 6);
reggrid=struct('startgrid', reggrid);


figrid0=startgrid;
figrid0.ABpts=linspace(0.03, 0.045, 6);
figrid0.MBpts=linspace(2.3, 3.1, 7);
figrid0.KREOpts=linspace(0.00, 0.03, 6);
figrid0.WIpts=[.01,.03,.05,.06,.08,.11,.13,.15,.19,.21,.23,.25];
figrid=struct('startgrid',figrid0);

noREOgrid=startgrid;
noREOgrid.KREOpts=linspace(0.00, 0.01, 6);
noREOgrid=struct('startgrid',noREOgrid);

asymgrid=startgrid;
asymgrid.ABpts=linspace(0.04, 0.08, 6);
asymgrid.MBpts=linspace(1.4, 2.6, 7);
asymgrid.KREOpts=linspace(0.00, 0.03, 6);
asymgrid.WIpts=[.03,.05,.06,.08,.11,.13,.15,.17,.19,.21,.23,.25,.3,.35];
asymgrid=struct('startgrid',asymgrid);

asymgrid2=startgrid;
asymgrid2.ABpts=linspace(0.027, 0.041, 6);
asymgrid2.MBpts=linspace(2.1, 2.9, 7);
asymgrid2.KREOpts=linspace(0, 0.03, 6);
asymgrid2.WIpts=[.03,.05,.06,.08,.11,.13,.15,.17,.19,.21,.23,.25,.3];
asymgrid2=struct('startgrid',asymgrid2);

onlygrid=startgrid;
onlygrid.ABpts=linspace(0.03, 0.043, 6);
onlygrid.MBpts=linspace(2.2, 2.9, 7);
onlygrid.KREOpts=linspace(0.00, 0.03, 6);
onlygrid.WIpts=[.03,.06,.1,.12,.15,.17,.2,.22,.24,.27,.3,.4,.6];
onlygrid=struct('startgrid',onlygrid);

onlygrid2=startgrid;
onlygrid2.ABpts=linspace(0.028, 0.042, 6);
onlygrid2.MBpts=linspace(2.3, 3.1, 7);
onlygrid2.KREOpts=linspace(0.00, 0.03, 6);
onlygrid2.WIpts=[-.02,.01,.03,.06,.1,.15,.17,.2,.22,.24,.27,.3];
onlygrid2=struct('startgrid',onlygrid2);

BGgrid=startgrid;
BGgrid.BGpts=[0,0.025,0.05,0.1];
BGgrid=struct('startgrid',BGgrid);

aggtaugrid=agggrid0;
aggtaugrid.BGpts=[0,0.025,0.05,0.08];
aggtaugrid.MBpts=linspace(2.0, 2.7, 7);
aggtaugrid=struct('startgrid', aggtaugrid);

liqgrid=startgrid;
liqgrid.ABpts=linspace(.035,.045,6);
liqgrid.MBpts=linspace(2.2,3.1,7);
liqgrid=struct('startgrid',liqgrid);

liqagggrid=agggrid0;
liqagggrid.ABpts=linspace(.032, .047, 7);
liqagggrid.MBpts=linspace(2.1, 3.0, 7);
liqagggrid=struct('startgrid', liqagggrid);

liqfigrid=figrid0;
liqfigrid.MBpts=linspace(2.3, 3.2, 7);
liqfigrid=struct('startgrid', liqfigrid);

liqreggrid=figrid0;
liqreggrid.WIpts=[-.06,-.03,-.01,.01,.03,.05,.06,.08,.11,.13,.15,.19,.20,.3,.35];
liqreggrid.MBpts=linspace(2.3, 3.2, 7);
liqreggrid=struct('startgrid', liqreggrid);

liqonlygrid=startgrid;
liqonlygrid.ABpts=linspace(0.035, 0.045, 6);
liqonlygrid.MBpts=linspace(2.3, 3.2, 7);
liqonlygrid.KREOpts=linspace(0.00, 0.03, 6);
liqonlygrid.WIpts=[-.06,-.03,-.01,.01,.03,.05,.06,.08,.11,.13,.15,.19,.20,.3,.35];
liqonlygrid=struct('startgrid',liqonlygrid);

liqonlygrid2=startgrid;
liqonlygrid2.ABpts=linspace(0.025, 0.038, 6);
liqonlygrid2.MBpts=linspace(1.9, 2.8, 7);
liqonlygrid2.KREOpts=linspace(0.00, 0.05, 6);
liqonlygrid2.WIpts=[-.03,-.01,.01,.03,.05,.06,.08,.11,.13,.15,.19,.20,.3,.35];
liqonlygrid2=struct('startgrid',liqonlygrid2);

liqasymgrid=startgrid;
liqasymgrid.ABpts=linspace(0.04, 0.06, 6);
liqasymgrid.MBpts=linspace(1.6, 2.8, 7);
liqasymgrid.KREOpts=linspace(0.00, 0.03, 6);
liqasymgrid.WIpts=[.03,.05,.06,.08,.11,.13,.15,.17,.19,.21,.23,.25,.3,.35];
liqasymgrid=struct('startgrid',liqasymgrid);

liqasymgrid2=startgrid;
liqasymgrid2.ABpts=linspace(0.03, 0.044, 6);
liqasymgrid2.MBpts=linspace(2.2, 3.1, 7);
liqasymgrid2.KREOpts=linspace(0.00, 0.05, 6);
liqasymgrid2.WIpts=[.03,.05,.06,.08,.11,.13,.15,.17,.19,.21,.23,.25,.3,.35];
liqasymgrid2=struct('startgrid',liqasymgrid2);



om25grad=linspace(0,0.25,5);
agggrad=linspace(0,1,5);

grids = struct('benchgrid',benchgrid, ...
    'agggrid',agggrid,...
    'sh3grid',sh3grid,...
    'reqgrid',reqgrid,...
    'reggrid',reggrid,...
    'figrid',figrid,...
    'asymgrid',asymgrid,...
    'asymgrid2',asymgrid2,...
    'onlygrid',onlygrid,...
    'onlygrid2',onlygrid2,...
    'noREOgrid',noREOgrid,...
    'BGgrid',BGgrid,...
    'aggtaugrid',aggtaugrid,...
    'liqgrid',liqgrid,...
    'liqagggrid',liqagggrid,...
    'liqfigrid',liqfigrid,...
    'liqreggrid',liqreggrid,...
    'liqasymgrid',liqasymgrid,...
    'liqasymgrid2',liqasymgrid2,...
    'liqonlygrid',liqonlygrid,...
    'liqonlygrid2',liqonlygrid2);


% Define experiments as modifications to base

expers_calib = {'bench',{};
                'om25',{'iota_om',0.25; 'iota_p',0; 'grid','figrid'};
                'agg',{'iota_om',0; 'iota_p',1; 'grid','agggrid'};
                'aggom25',{'iota_om',0.25; 'iota_p',1; 'grid','reggrid'};
                'aggom25ABonly',{'iota_om',0.25; 'iota_p',1; 'indexMB',0; 'grid','onlygrid'};
                'aggom25MBonly',{'iota_om',0.25; 'iota_p',1; 'indexAB',0; 'grid','onlygrid2'};
                'aggom25max100',{'pgr_max',1.0; 'iota_om',0.25; 'iota_p',1; 'grid','asymgrid'};
                'aggom25ABonlymax100',{'pgr_max',1.0;'iota_om',0.25; 'iota_p',1; 'indexMB',0;  'grid','asymgrid2'};
                'aggomTail25',{'iota_om',0.25; 'iota_p',1; 'indexTail',1; 'pgr_max',0.9; 'grid','liqasymgrid'};
                'aggom25s1',{'iota_om',om25grad(2); 'iota_p',agggrad(2); 'grid','benchgrid'};
                'aggom25s2',{'iota_om',om25grad(3); 'iota_p',agggrad(3); 'grid','benchgrid'};
                'aggom25s3',{'iota_om',om25grad(4); 'iota_p',agggrad(4); 'grid','reggrid'};
                'aggom25o1',{'iota_om',om25grad(2); 'iota_p',1; 'grid','agggrid'};
                'aggom25o2',{'iota_om',om25grad(3); 'iota_p',1; 'grid','reggrid'};
                'aggom25o3',{'iota_om',om25grad(4); 'iota_p',1; 'grid','reggrid'};
                'aggom25a1',{'iota_om',.25; 'iota_p',agggrad(2); 'grid','figrid'};
                'aggom25a2',{'iota_om',.25; 'iota_p',agggrad(3); 'grid','reggrid'};
                'aggom25a3',{'iota_om',.25; 'iota_p',agggrad(4); 'grid','reggrid'};
                'nosav',{'psiS0',1e5};
		        'nosavagg',{'psiS0',1e5; 'iota_p',1; 'grid','agggrid'};
                'nosavom25',{'psiS0',1e5; 'iota_om',0.25; 'grid','figrid'};
                'nosavaggom25',{'psiS0',1e5; 'iota_p',1; 'iota_om',0.25; 'grid','reggrid'}                
                'shI3aggom25',{'chiI',0.03;'shareI',.0385; 'shareHI',.0385; 'iota_p',1; 'iota_om',0.25; 'grid','sh3grid'};
                'shI4aggom25',{'chiI',0.04;'shareI',.051; 'shareHI',.051; 'iota_p',1; 'iota_om',0.25; 'grid','reggrid'}                
                'shI6aggom25',{'chiI',0.06;'shareI',.076; 'shareHI',.076; 'iota_p',1; 'iota_om',0.25; 'grid','figrid'};                            
                'shI10aggom25',{'chiI',0.1;'shareI',.126; 'shareHI',.126; 'iota_p',1; 'iota_om',0.25; 'grid','figrid'}                
                'taxL080',{'tauL',0.8; 'grid','BGgrid'};                
                'aggtauL080',{'tauL',0.8; 'iota_om',0; 'iota_p',1; 'grid','aggtaugrid'};
                'aggcapreq90',{'iota_om',0; 'iota_p',1; 'phiI',0.9; 'grid','reqgrid'};
                'aggcapreq85',{'iota_om',0; 'iota_p',1; 'phiI',0.85; 'grid','reqgrid'};               
                'liqhi',{'theta',0.18;'etaB',.05;'liq_thresh',0.9;'sigom_lo',0.16;'sigom_hi',0.21;'xilo',.13;'xihi',.19; 'indexMB',0; 'indexAB',0; 'iota_om',0; 'grid','liqgrid'};
		        'liqhiagg',{'theta',0.18;'etaB',.05;'liq_thresh',0.9;'sigom_lo',0.16;'sigom_hi',0.21;'xilo',.13;'xihi',.19; 'iota_p',1; 'grid','liqagggrid'};
                'liqhiom25',{'theta',0.18;'etaB',.05;'liq_thresh',0.9;'sigom_lo',0.16;'sigom_hi',0.21;'xilo',.13;'xihi',.19; 'iota_om',0.25; 'grid','liqfigrid'};
                'liqhiaggom25',{'theta',0.18;'etaB',.05;'liq_thresh',0.9;'sigom_lo',0.16;'sigom_hi',0.21;'xilo',.13;'xihi',.19; 'iota_p',1; 'iota_om',0.25; 'grid','liqreggrid'};                
                'liqhiaggom25ABonly',{'theta',0.18;'etaB',.05;'liq_thresh',0.9;'sigom_lo',0.16;'sigom_hi',0.21;'xilo',.13;'xihi',.19; 'iota_p',1; 'iota_om',0.25; 'indexMB',0; 'grid','liqgrid'};                
                'liqhiaggom25MBonly',{'theta',0.18;'etaB',.05;'liq_thresh',0.9;'sigom_lo',0.16;'sigom_hi',0.21;'xilo',.13;'xihi',.19; 'iota_p',1; 'iota_om',0.25; 'indexAB',0; 'grid','liqonlygrid'};                
                'liqhiaggom25max100',{'theta',0.18;'etaB',.05;'liq_thresh',0.9;'sigom_lo',0.16;'sigom_hi',0.21;'xilo',.13;'xihi',.19; 'iota_p',1; 'iota_om',0.25; 'pgr_max',1.0; 'grid','liqasymgrid'};                
                'liqhiaggom25ABonlymax100',{'theta',0.18;'etaB',.05;'liq_thresh',0.9;'sigom_lo',0.16;'sigom_hi',0.21;'xilo',.13;'xihi',.19; 'iota_p',1;  'iota_om',0.25; 'pgr_max',1.0; 'indexMB',0; 'grid','liqasymgrid2'};                
                'liqhiaggomTail25',{'theta',0.18;'etaB',.05;'liq_thresh',0.9;'sigom_lo',0.16;'sigom_hi',0.21;'xilo',.13;'xihi',.19; 'iota_p',1; 'iota_om',0.25; 'indexTail',1; 'pgr_max',0.9; 'grid','liqasymgrid'} ;               
                };

% Second, vary other parameters
expers_other = { };
expers = [expers_calib; expers_other];


%% Create each experiment definition (fully dynamic below this line)
N_EXP=size(expers,1);
N_GRIDS = length(fieldnames(benchgrid));

% Write list of defined experiments to file
fid = fopen([mfilename,'.txt'],'w');
for i=1:N_EXP
   fprintf(fid,expers{i,1});
   fprintf(fid,'\n');
end
fclose(fid);

% Create experiment table
baserow = {'test','ini0',{startgrid.ABpts},{startgrid.MBpts},{startgrid.KREOpts},{startgrid.WIpts},{startgrid.BGpts}};
basetable = cell2table(baserow);          
basetable.Properties.VariableNames = {'exper','grid','ABpts','MBpts','KREOpts','WIpts','BGpts'};
basetable = [basetable, struct2table(params,'AsArray',true)];

expertable = repmat(basetable,N_EXP*N_GRIDS,1);

benchgrid_cell=struct2cell(benchgrid);
for i=1:N_GRIDS
    fnames = fieldnames(benchgrid_cell{i});
    gridvectors=struct2cell(benchgrid_cell{i})';
    [~,~,order]=intersect(expertable.Properties.VariableNames(3:7),fnames, ...
            'stable');
    expertable{N_EXP*(i-1)+1 : N_EXP*i,3:7} = ...
        repmat( gridvectors(order), N_EXP, 1);
    expertable(N_EXP*(i-1)+1 : N_EXP*i,2) = ...
        {benchgrid_cell{i}.label};   
end

gridOrder = @(grid_as_cell,order)grid_as_cell(order);
for i=1:N_EXP
   expertable.exper([i,N_EXP+i]) = expers(i,1);
   changes = expers{i,2};
   for j=1:size(changes,1)
      varname=changes{j,1};
      if strcmp(varname,'grid')
          newgrid = struct2cell(grids.(changes{j,2}));
          fnames = fieldnames(newgrid{1});
          [~,~,order]=intersect(expertable.Properties.VariableNames(3:7),fnames, ...
            'stable');
        for k=1:N_GRIDS
          expertable{i+N_EXP*(k-1),3:7} = ...
              gridOrder(struct2cell(newgrid{k}),order)';
        end
      else
          if ischar(changes{j,2})
            expertable.(varname)([i,N_EXP+i],:)=repmat(changes(j,2),N_GRIDS,1);
          else
            expertable.(varname)([i,N_EXP+i],:)=repmat(changes{j,2},N_GRIDS,1);    
          end
      end
   end
end

% Name each row
row_id = cell(2*N_EXP,1);
idx_suffix = find(~strcmp(expertable.grid,''));
idx_nosuffix = find(strcmp(expertable.grid,''));
row_id(idx_suffix) = arrayfun(@(idx){[expertable.exper{idx},'_',expertable.grid{idx}]},idx_suffix,'UniformOutput',false);
row_id(idx_nosuffix) = arrayfun(@(idx)expertable.exper(idx),idx_nosuffix,'UniformOutput',false);
expertable.Properties.RowNames = [row_id{:}]';
%expertable.Properties.RowNames(1) = {'xpbase'};

% Create a struct of structs
allexpers=cell2struct(expertable.Properties.RowNames,expertable.Properties.RowNames);
for i=1:size(expertable,1)
    s = table2struct(expertable(i,3:7));
    s.params = table2struct(expertable(i,8:end));
    name = expertable.Properties.RowNames(i);
    name = name{:};
    allexpers.(name) = s;
end

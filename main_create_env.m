if ~exist('clear_flag', 'var'), clear_flag = 1; end

if usejava('desktop') && clear_flag
   clear;
end

% ===========================================
% program control
% ===========================================

% name of file with experiment definitions
if ~exist('experdef_file','var')
    experdef_file='experdef_20191112';
end
% name of experiment
if ~exist('expername','var')
     expername='bench_ini0';
%    expername='liqhi_ini0';
end

% possible values 'no_guess', 'guess'
if ~exist('guess_mode','var')
    guess_mode='no_guess';
end
% path to file with initial guess; not used if guess_mode='no_guess'
if ~exist('guess_path','var')
    guess_path='res_20190828_shI4aggom25_i150c100.mat';
end

% path and file name for output file
outfname=['env_',expername,'.mat'];

% approximation mode
if ~exist('approxMode','var')
    approxMode='linear';
end

% Use analytic Jacobian (requires Symbolic Math Toolbox)
% set to false if not installed
% not creating Jacobian is slow, does not work well with Matlab 2019a or
% higher
if ~exist('useJacobian','var')
    useJacobian=true;
end

% only compute steady state (without writing initial model object)
if ~exist('ststonly', 'var'), ststonly=0; end

% ===================================================
% non-stochastic steady state
% ===================================================

% load experiment definitions
run(experdef_file);
disp(experdef_file);
expdef=allexpers.(expername);


% Update parameters
params = expdef.params;
params.chiS=1-params.chiB-params.chiI;
params.shareS=1-params.shareB-params.shareI;
params.shareHS=1-params.shareHB-params.shareHI;
expdef.params = params;
% End of paramter update
modelclassname='SAMIntermediaryModel';

ststfun=str2func([modelclassname,'.compStSt']);
instfun=str2func(modelclassname);
guessfun=str2func([modelclassname,'.assignGuess']);
jacfun=str2func([modelclassname,'.constructJacobian']);

                        
% compute steady state values
% do this once for each state of sigma_omega
options=optimset('Display','iter','TolX',1e-10,'TolFun',1e-10,'MaxIter',100);


% No indexation version
gvec = {
    [ 2.6920   -8.9381   -1.0239   -2.1973   -1.0046   -0.6184    0.3231    1.2126   -0.8188    1.4397   -2.6597],...
    [ 2.6920   -8.9381   -1.0239   -2.1973   -1.0046   -0.6184    0.3231    1.2126   -0.8188    1.4397   -2.6597]
};


stv=cell(2,1);
for i=1:2
	fh_compStSt=@(x)ststfun(x,expdef.params,i,0);
	[solvec,~,exfl]=fsolve(fh_compStSt,gvec{i},options);
	if exfl<1
		disp('!! Problem computing steady state');
	end
    [~,stvtmp]=ststfun(solvec,expdef.params,i,1);
	stv{i}=stvtmp;
   
end

default_rate = 0.75 * (1.0 - stv{1}.statsout.ZN ) + 0.25 * (1.0 - stv{2}.statsout.ZN);
fprintf('Total default rate: %4.3f\n', default_rate);
fprintf('Fraction of liquidity defaults (boom): %4.3f\n', stv{1}.statsout.frac_liq);
fprintf('Fraction of liquidity defaults (bust): %4.3f\n', stv{2}.statsout.frac_liq);
fprintf('Strategic default penalty: %4.3f\n', stv{1}.statsout.strat_pen);

expdef.params.KBbar=stv{1}.statsout.KBbar;
expdef.params.HI=stv{1}.statsout.HI;
expdef.params.HS=stv{1}.statsout.HS;

if ststonly
    return;
end


% ===================================================
% Parameters of stochastic model
% ===================================================

% for mean-reverting AR(1) (log)TFP, denoted x here
N_y=5;
sig_y=expdef.params.sig_y;
mu_y=expdef.params.mu_y;
mu_Y=expdef.params.mu_Y;
rho_y=expdef.params.rho_y;
Skew_y=0;
[Yprob,y] = model.DSGEModel.rouwen(rho_y,mu_y,sig_y,Skew_y,N_y);
Y=exp(y); % exponentiated to get TFP from log(TFP)
Yprob=Yprob';

% for shock to std.dev. of borrowers house values
N_om=2;
sig2_om=log(1+[expdef.params.sigom_lo,expdef.params.sigom_hi].^2);
expdef.params.sig2_om=sig2_om;
Om=sig2_om';
xivec=[expdef.params.xihi,expdef.params.xilo];
expdef.params.xi=xivec';


Omprob_recess=reshape(expdef.params.Omprob_recess,2,2);
Omprob_boom=reshape(expdef.params.Omprob_boom,2,2);

% total Markov transition for exogenous states
% first find recession states

threshold=mu_y;
recess=find(Y<1+threshold,1,'last');
boom=N_y-recess;
omcomb=grid.StateSpaceGrid.makeCombinations([N_y,N_om]);
mpts_perm=[ Y(omcomb(:,1)), Om(omcomb(:,2)), expdef.params.xi(omcomb(:,2)) ];

% transition matrix
rectrans=kron(Yprob(1:recess,:),Omprob_recess);
boomtrans=kron(Yprob(recess+1:end,:),Omprob_boom);
mtrans=[ rectrans; boomtrans ];
exnpt=size(mpts_perm,1);
mpts_all=[mpts_perm, (1:exnpt)'];
     
% simulate exog process and produce histogram
Nsim=100000;
simst=model.DSGEModel.hitm_s(mtrans',rand(Nsim,1));
simfrac=histc(simst,1:exnpt)/Nsim;
disp('-----------------------------------------');
disp('Unconditional prob of each state: ');
disp(num2str([(1:exnpt)',mpts_perm,simfrac]));
disp(['"Recession" states: ', num2str(sum(simfrac(1:N_om*recess)))]); 
omhi_ind=N_om*(1:N_y);
disp(['High uncertainty states: ', num2str(sum(simfrac(omhi_ind)))]);
% average length of spell
forecl=ismember(simst,omhi_ind);
forestend=forecl(2:end)-forecl(1:end-1);
firstend=find(forestend==-1,1,'first');
laststart=find(forestend==1,1,'last');
forestend=forestend(firstend+1:laststart-1);
forecep=sum(forestend==1);
disp(['Uncertainty episodes (%): ',num2str(forecep/Nsim)]);
forest=find(forestend==1);
foreend=find(forestend==-1);
forelen=foreend-forest;
disp(['Average length: ',num2str(mean(forelen))]);


% variable lists
% exogenous state variable names
exogenv=struct;
exogenv.exnames={'Y','Om','xi'};
exogenv.exidx=struct('ZA',1,'Om',2,'xi',3);
exogenv.exnpt=exnpt;
exogenv.pts_perm=mpts_perm;
exogenv.pts_all=mpts_all;
exogenv.mtrans=mtrans;
exogenv.normparam=[mu_y, sig_y^2];

          
% ===========================================
% create model object
% ===========================================

% endogenous state variables
endogenv=struct;

if expdef.params.consolidatedLender
    endogenv.ennames={'AB','MB','KREO'};
    % initial definition of grids based on steady state values of state vars
    unigrids={1:exogenv.exnpt,expdef.ABpts,expdef.MBpts,expdef.KREOpts};    
else
    endogenv.ennames={'AB','MB','KREO','WI','BG'};
    % initial definition of grids based on steady state values of state vars
    unigrids={1:exogenv.exnpt,expdef.ABpts,expdef.MBpts,expdef.KREOpts,expdef.WIpts,expdef.BGpts};
end

% assign values for initial guess
[solguessvec,Vguessvec,Vnames]=guessfun(stv{1});


% save name lists
endogenv.solnames=fieldnames(stv{1}.Sol);
endogenv.solbase=solguessvec;
endogenv.addnames=fieldnames(stv{1}.Add);
endogenv.Vnames=Vnames;
endogenv.condnames = {'expRP_A', 'expRP_M', 'expRP_T', 'CovP_A', 'Del_BC_K', 'Del_BC_M'};                        
basegrid=grid.TensorGrid(unigrids);
            

if isequal(guess_mode,'no_guess')
    inigrid=basegrid;
    solmatguess=repmat(solguessvec',[basegrid.Npt,1]);
    Vmatguess=repmat(Vguessvec',[basegrid.Npt,1]);
    transguess=kron(basegrid.Pointmat(:,2:end), ones(1,exogenv.exnpt));
else
    %inigrid=grid.TensorGrid(basegrid.StateBounds,[exnpt,20,20,15]);
    inigrid=basegrid;
    guessobj=load(guess_path,'mobj');
    solmatguess=guessobj.mobj.evaluatePol(inigrid.Pointmat)';
    Vmatguess=guessobj.mobj.evaluateVal(inigrid.Pointmat)';
    transguess=guessobj.mobj.evaluateTrans(inigrid.Pointmat)';
    
end                  
            
% build approximating function
% create guess for solution variable functions
Pf=grid.LinearInterpFunction(inigrid,solmatguess);
% create guess for next-period functions
Vf=grid.LinearInterpFunction(inigrid,Vmatguess);
% and for state transition function
Tf=grid.LinearInterpFunction(inigrid,transguess);
% and for temporary value function (for risk-taker
% bankruptcy)
Tempf=grid.LinearInterpFunction(inigrid,zeros(inigrid.Npt,1));


if useJacobian
   J=jacfun(expdef.params); 
else
   J=[];
end
            
% create new object
mobj=instfun(expdef.params,endogenv,exogenv,Vf,Pf,Tf,Tempf,basegrid,J);

disp('-------------------------------------');
disp('Bounds of state space:');
disp(num2str(mobj.Vfct.SSGrid.StateBounds));

% ===========================================
% save initial object
% ===========================================

save(outfname,'stv','mobj');              







if ~exist('clear_flag', 'var'), clear_flag = 1; end

if usejava('desktop') && clear_flag
   clear;
end
close all;

% ===========================================
% program control
% ===========================================

open_parpool; % HERE

% path to file with experiment definition
if ~exist('exper_path','var')
    exper_path='env_bench_ini0.mat';
end
if ~exist('maxit','var')
    maxit=150;
end
if ~exist('tol_avg','var')
    % mean convergence
    tol_avg=1e-5;
end


if ~exist('maxit_VF','var')
    % extra VF convergence iterations?
    maxit_VF = 100;
end

% print mode
printmode=0;

% dampening
damp=0.5;


% ===========================================
% load environment structure
% ===========================================

load(exper_path);

% ===========================================
% outer convergence loop
% ===========================================

% policy iteration
[mobj,failedPoints,dist,distT]=mobj.polIter(maxit,1,printmode,damp,tol_avg);

if maxit_VF>0
    [mobj,distVF,mean_distVF]=mobj.convergeVF(maxit_VF,0,tol_avg);
end

% make date string
curr_date=clock;
date_str='';
for i=1:5
    date_str=[date_str,'_',num2str(curr_date(i))];
end

if ~exist('outname','var')
    outname=['res',date_str,'.mat'];
else
    outname=[outname,date_str,'.mat'];    
end

save(outname,'mobj','stv','failedPoints','dist','distT');



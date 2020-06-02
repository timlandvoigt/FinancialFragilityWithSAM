function transition_value()

%--------------------------------------------------------------------------
% INITIAL SETUP
%--------------------------------------------------------------------------

if ~exist('clear_flag', 'var'), clear_flag = 1; end

if usejava('desktop') && clear_flag
    clear;
else
    ver;
end
close all;

no_par_processes = 4;
open_parpool;

respath='./';
outpath='./Results/';
suffix = '';
recompute = 0;

% Baseline file (simulations will come from here)
if ~exist('resfile1','var')
    resfile1 = ['res_20180628_rebmaint', suffix];
end

% Transition files (value functions will come from here)
if ~exist('resfile2_list','var')
    resfile2_list = { ...
        'res_20180628_rebagg', ...
        'res_20180628_rebom25', ...
        'res_20180628_rebaggom25', ...
        'res_20180628_rebaggom25max100', ...
        };
end

%         'res_20180628_bench', ...
%         'res_20180628_agg', ...
%         'res_20180628_om25', ...
%         'res_20180628_aggom25', ...
%         'res_20180628_aggom25ABonly', ...
%         'res_20180628_aggom25ABonlymax100', ...
%         'res_20180628_aggom25max100', ...
%         'res_20180628_aggom25MBonly', ...

print_output = 1;

%--------------------------------------------------------------------------
% GET SIMULATED POINTS FROM MODEL 1
%--------------------------------------------------------------------------

% Load baseline economy
load([respath,resfile1,'.mat']);

% Initial Economy Config
varlist={'simseries','statevec','indexmap','varnames'};
load(['sim_',resfile1],varlist{:});

% set starting point
statevec=statevec(2:end);

% number of periods and burn-in
N_runs = 10000;
N_sim = 10000;

% Get index for endogenous states
if mobj.Params.consolidatedLender
%     envarind=4:6;
    envarind=4:7;
else
%     envarind=4:7;
    envarind=4:8; % Adding government debt
end

% Evaluate value function simulated states
startptmat = [statevec, simseries(:, envarind)];
startptmat = startptmat(1 : N_sim, :);

values1 = simseries(1 : size(startptmat, 1), :);

steady_names1 = varnames;
steady_mean1 = mean(simseries);
steady_std1 = std(simseries);

% values1 = mobj.evaluateVal(startptmat);
% soln1 = mobj.evaluatePol(startptmat);
% 
% % Compute moments
% mean_vals1 = mean(values1, 2);
% std_vals1 = std(values1, 0, 2);

% simseries2 = zeros(size(simseries));
% for jj = 1 : N_runs
%     
% end

% First one to get sizes

[test, ~, ~, ~] = compute_trans_full(mobj, startptmat(1, :));
ix_good = abs(test - simseries(1, :)) < 1e-8;

tex_dir = ['./tex/values/', resfile1, '/'];
if ~exist(tex_dir, 'dir')
    mkdir(tex_dir);
end

% HERE
% Steady state version
for ii = 1 : length(varnames)
    if ix_good(ii)
        name_i = varnames{ii};

        outfile = sprintf('%s%s_steady_level.tex', tex_dir, name_i);
        fid = fopen(outfile, 'w');
        fprintf(fid, '%4.3f', steady_mean1(ii));
        fclose(fid);

        outfile = sprintf('%s%s_steady_pct_level.tex', tex_dir, name_i);
        fid = fopen(outfile, 'w');
        fprintf(fid, '%3.2f', 100.0 * steady_mean1(ii));
        fclose(fid);
    end
end

if print_output
    [names_sorted, names_ix] = sort(varnames);
    fprintf('\n\nMODEL 1 STEADY:\n\n');
    for ii = 1 : length(varnames)
        fprintf('%s: %4.3g\n', names_sorted{ii}, steady_mean1(names_ix(ii)));
    end
end

% values = values1;
% save(['trans_', resfile1], 'values', 'varnames');

%--------------------------------------------------------------------------
% NOW EVALUATE ON MODEL 2
%--------------------------------------------------------------------------

for j_file = 1 : length(resfile2_list)
    
    tic
    
    % Load comparison model
    resfile2 = [resfile2_list{j_file}, suffix];
    fprintf('Comparison model: %s\n', resfile2);
    
    tex_dir = ['./tex/values/', resfile2, '/'];
    if ~exist(tex_dir, 'dir')
        mkdir(tex_dir);
    end
        
    load([respath,resfile2,'.mat']);
    varlist={'simseries','statevec','indexmap','varnames'};
    load(['sim_',resfile2],varlist{:});
    
    % Steady state moments
    steady_names2 = varnames;
    steady_mean2 = mean(simseries);
    steady_std2 = std(simseries);
    
%     values2 = mobj.evaluateVal(startptmat);
    
    trans_file = ['trans_', resfile2];
    if recompute || ~exist([trans_file, '.mat'], 'file')
        [values2, names2, mean2, std2] = compute_trans_full(mobj, startptmat);
        values = values2;
        save(trans_file, 'values', 'varnames');
    else
        [values2, names2, mean2, std2] = load_trans_file(trans_file);
    end
    
    % Compute moments
%     mean_vals2 = mean(values2, 2);
%     std_vals2 = std(values2, 0, 2);
    
    % Compute differences
%     diff_vals = values2 - values1;
    mean_vals_diff = mean(values2, 1) - steady_mean1;
    mean_pct_diff = 100.0 * (mean_vals_diff ./ steady_mean1);
    
    steady_mean_vals_diff = steady_mean2 - steady_mean1;
    steady_mean_pct_diff = 100.0 * steady_mean_vals_diff ./ steady_mean1;
    
%     std_vals_diff = std(diff_vals, 0, 2);
%     std_pct_diff = std_vals_diff ./ mean_vals1;
    
%     if print_output
%         fprintf('\n\nDIFFERENCE:\n\n');
%         for ii = 1 : mobj.NV
%             fprintf('%s: %4.3g%%\n', mobj.V_names{ii}, mean_pct_diff(ii));
%         end
%     end
    toc    

    if print_output
        [names_sorted, names_ix] = sort(names2);
        
        fprintf('\n\nMODEL 1 STEADY:\n\n');
        for ii = 1 : length(names2)
            fprintf('%s: %3.2f\n', names_sorted{ii}, steady_mean1(names_ix(ii)));
        end
        
        fprintf('\n\nMODEL 2 TRANSITION:\n\n');
        for ii = 1 : length(names2)
            fprintf('%s: %3.2f\n', names_sorted{ii}, mean2(names_ix(ii)));
        end
        
        fprintf('\n\nMODEL 2 STEADY:\n\n');
        for ii = 1 : length(names2)
            fprintf('%s: %3.2f\n', names_sorted{ii}, steady_mean2(names_ix(ii)));
        end
        
        fprintf('\n\nMODEL 2 STEADY PCT DIFF:\n\n');
        for ii = 1 : length(names2)
            fprintf('%s: %3.2f\n', names_sorted{ii}, steady_mean_pct_diff(names_ix(ii)));
        end
    end
    
    % Save output files
    for ii = 1 : length(varnames)
        if ix_good(ii)
            name_i = varnames{ii};

            outfile = sprintf('%s%s_level.tex', tex_dir, name_i);
            fid = fopen(outfile, 'w');
            fprintf(fid, '%4.3f', mean2(ii));
            fclose(fid);
            
            outfile = sprintf('%s%s_level_pct.tex', tex_dir, name_i);
            fid = fopen(outfile, 'w');
            fprintf(fid, '%3.2f', 100.0 * mean2(ii));
            fclose(fid);

            outfile = sprintf('%s%s_abs_diff.tex', tex_dir, name_i);
            fid = fopen(outfile, 'w');
            fprintf(fid, '%+4.3f', mean_vals_diff(ii));
            fclose(fid);
            
            outfile = sprintf('%s%s_ppt_diff.tex', tex_dir, name_i);
            fid = fopen(outfile, 'w');
            fprintf(fid, '%+3.2f', 100.0 * mean_vals_diff(ii));
            fclose(fid);

            outfile = sprintf('%s%s_pct_diff.tex', tex_dir, name_i);
            fid = fopen(outfile, 'w');
            fprintf(fid, '%+3.2f', mean_pct_diff(ii));
            fclose(fid);
        end
    end
    
    % Steady state version
    for ii = 1 : length(varnames)
        if ix_good(ii)
            name_i = varnames{ii};

            outfile = sprintf('%s%s_steady_level.tex', tex_dir, name_i);
            fid = fopen(outfile, 'w');
            fprintf(fid, '%4.3f', steady_mean2(ii));
            fclose(fid);
            
            outfile = sprintf('%s%s_steady_pct_level.tex', tex_dir, name_i);
            fid = fopen(outfile, 'w');
            fprintf(fid, '%3.2f', 100.0 * steady_mean2(ii));
            fclose(fid);

            outfile = sprintf('%s%s_steady_abs_diff.tex', tex_dir, name_i);
            fid = fopen(outfile, 'w');
            fprintf(fid, '%+4.3f', steady_mean_vals_diff(ii));
            fclose(fid);
            
            outfile = sprintf('%s%s_steady_ppt_diff.tex', tex_dir, name_i);
            fid = fopen(outfile, 'w');
            fprintf(fid, '%+3.2f', 100.0 * steady_mean_vals_diff(ii));
            fclose(fid);

            outfile = sprintf('%s%s_steady_pct_diff.tex', tex_dir, name_i);
            fid = fopen(outfile, 'w');
            fprintf(fid, '%+3.2f', steady_mean_pct_diff(ii));
            fclose(fid);
        end
    end
    
end

end

function [transseries_full, transnames_full, mean_trans, std_trans] = compute_trans_full(mobj, startptmat)

    N_runs = size(startptmat, 1);

    [temp, transnames, ~, ~, ~] = mobj.simulate(1, 0, startptmat(1, :), 0, int64(1));
    N_trans = length(temp);
    transseries = zeros(N_runs, N_trans);
    transseries(1, :) = temp;
    parfor jj = 2 : N_runs
        [transseries(jj, :), ~, ~, ~, ~] = mobj.simulate(1, 0, startptmat(jj, :), 0, int64(1));
    end

    stacked_series_j = [NaN * ones(1, N_trans); transseries(1, :)];
    [temp, transnames_full] = mobj.computeSimulationMoments(stacked_series_j, transnames);
    transseries_full = zeros(N_runs, length(temp));
    transseries_full(1, :) = temp;
    parfor jj = 2 : N_runs
        stacked_series_j = [NaN * ones(1, N_trans); transseries(1, :)];
        [transseries_full(jj, :), ~] = mobj.computeSimulationMoments(stacked_series_j, transnames);
    end
    
    mean_trans = mean(transseries_full);
    std_trans = std(transseries_full);

end

function [values, varnames, means, stds] = load_trans_file(trans_file)

load(trans_file);
means = mean(values);
stds = std(values);

end

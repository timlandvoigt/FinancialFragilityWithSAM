results_files = dir('sim_res_20190112_*.mat');

pzns_res = 'res_20191112_bench';
outfolder = './Results/';

bench.pzns=load(['pzns_',pzns_res]);
bench.res=load(pzns_res);

values=struct2cell(bench.res.mobj.Params);
paramout=cell2table(values,'RowNames',fieldnames(bench.res.mobj.Params)); 

EXOG_STATE=5;
state_range=4:8;
state = [EXOG_STATE, mean(bench.pzns.simseries(:,state_range))];
val_bench = bench.res.mobj.evaluateVal(state);
VB_bench = val_bench(4);
VI_bench = val_bench(5);
VS_bench = val_bench(6);
pzns_bench = bench.pzns.Pricefct.evaluateAt(state);
qCB_bench = pzns_bench(2);
qCI_bench = pzns_bench(3);
qCS_bench = pzns_bench(4);

bench = [];

EVs=zeros(size(results_files,1),7);
EVs=array2table(EVs);
EVs.Properties.VariableNames={'pctdiff_VB','pctdiff_VI','pctdiff_VS','deltaWB','deltaWI','deltaWS','Total'};
for i=1:size(results_files,1)
   fname = results_files(i).name;
   parts = strread(fname,'%s','delimiter','_');
   alt.sim=load(fname,'simseries','errmat','indexmap');
   alt.res=load(fname(5:end),'mobj');
   % For Each alternative economy, compute VX_bench / VX_alternative
   % diff_VX = how much higher agent's welfare is in alternative vs. bench
   % diff_VX = by how many percent do you need to increase consumption in 
   % every state of the world in benchmark economy to make her as well 
   % of as she was in the alternative
   alt_state =  [EXOG_STATE, mean(alt.sim.simseries(:,state_range))];
   val_alt = alt.res.mobj.evaluateVal(alt_state);
   VB_alt = val_alt(4);
   VI_alt = val_alt(5);
   VS_alt = val_alt(6);
   diff_VB = VB_alt / VB_bench - 1;
   diff_VI = VI_alt / VI_bench - 1;
   diff_VS = VS_alt / VS_bench - 1;
   
   EVs{i,1} = diff_VB;
   EVs{i,2} = diff_VI;
   EVs{i,3} = diff_VS;
   EVs{i,4} = diff_VB * qCB_bench;
   EVs{i,5} = diff_VI * qCI_bench;
   EVs{i,6} = diff_VS * qCS_bench;
   EVs.Properties.RowNames(i) = parts(4);
   
end
EVs{:,7} = EVs{:,4} + EVs{:,5} + EVs{:,6};

for j=1:4
    writetable(EVs,[outfolder,'statsexog_res_20191112_welfare.xls'],'Sheet',j,'WriteRowNames',true);
end
writetable(paramout,[outfolder,'statsexog_res_20191112_welfare.xls'],'Sheet','params','WriteRowNames',true);
writetable(EVs,[outfolder,'welfare_',pzns_res,'.xls'],'WriteRowNames',true);


% parallel processes? set this to zero for just one process
if ~exist('no_par_processes','var')
    no_par_processes=16; %10
end

if usejava('desktop')
    disp('Local Job')
    mycluster = 'local';
else
    
    if contains(getenv('QUEUE'), 'aws-')
        disp('AWS Job')
        mycluster = 'local';
        %myfolderbase = strcat(getenv('TMP'), '/');
    else
        disp('HPCC Job')
        mycluster = 'WhartonHPC';
        %myfolderbase = '~/matlabtmp/';
    end
    
end

mypool = parcluster(mycluster);

disp(['PPN: ',num2str(no_par_processes)]);

cp=gcp('nocreate');
if ~isempty(cp)
    if  cp.NumWorkers~=no_par_processes
        delete(cp);
        if no_par_processes>0
            parpool(mypool, no_par_processes);
        end
    end
else
    if no_par_processes>0
        parpool(mypool, no_par_processes);
    end
end

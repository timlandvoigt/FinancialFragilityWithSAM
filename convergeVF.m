clear;
close all;

clear_flag=0;
% name of file with experiment definitions
experdef_file='experdef_20190828';
% name of experiment
expername='aggom25MBonly_ini0';
% possible values 'no_guess', 'guess'
guess_mode='guess';
% path to file with initial guess; not used if guess_mode='no_guess'
guess_path='res_20190828_aggom25MBonly_i150';
% approximation mode
useJacobian=false;
% only compute steady state (without writing initial model object)
ststonly=0;


main_create_env;

stateSpace=mobj.Pfct.SSGrid.Pointmat;
npt=mobj.Pfct.SSGrid.Npt;
NDIM=mobj.Vfct.SSGrid.Ndim;
exnpt=size(mobj.Exogenv.pts_perm,1);

open_parpool;

vfindex=4:6;
tolvf=1e-6;

iter=0;
maxiter=500;

transmat=mobj.evaluateTrans(stateSpace)';

% Evaluate value functions at transition points
transpts=reshape([repmat(1:exnpt,npt,1),transmat],npt*exnpt,NDIM);
% index matching for transitions
refidx=kron((1:npt)',ones(1,exnpt));
refidx=refidx(:);
polmat=mobj.evaluatePol(stateSpace)';

while true
    
    iter = iter+1;
    Vtrans = mobj.evaluateVal(transpts)';
    oldV = mobj.evaluateVal(stateSpace)';
    newV = zeros(size(oldV));
    
    parfor p=1:npt
        
        thispt=stateSpace(p,:);
        thispol=polmat(p,:);

        thistrans=transmat(p,:)';
        vtrans=Vtrans(refidx==p,:);
        thisvals=oldV(p,:);        
                
        [nextst,outstr]=mobj.calcStateTransition(thispt,thispol,0,thistrans,thisvals);        
        p0=exp(thispol(4));
        thisexp = mobj.computeExpectations(thispt(1),nextst, p0, vtrans );

        %calcEquations(obj,exst,nextst,solvec,instr,mode,varargin)
        [~,~,thisV]=mobj.calcEquations(thispt(1),nextst,thispol,outstr,1,thisexp);
        
        newV(p,:)=thisV{1};
                
    end
        
    val = oldV(:,vfindex);
    nextval = newV(:,vfindex);
    [dist,wh]=max(abs(val(:) - nextval(:)));
    [mean_dist,col]=max(abs(mean(val-nextval)));
    [wh_1,wh_2]=ind2sub(size(nextval),wh);
    disp(['-- Iteration: ',num2str(iter),', max distance: ',num2str(dist),' in ',char(mobj.V_names(vfindex(1)-1+wh_2)), ...
        ' at point ',num2str(wh_1),': ',num2str(stateSpace(wh_1,:))]);
    disp(['-- Mean distance: ',num2str(mean_dist),' in ',char(mobj.V_names(vfindex(1)-1+col))]);
    disp('----------------------------------------------');
    if mean_dist < tolvf
        disp('Converged');
        break;
    end
    if iter >= maxiter
        disp('Max.iter. exceeded');
        break;
    end
    
    mobj=mobj.updateVfct(newV);

    
end





clear;
close all;
% Price an asset in zero net supply using already computed consumption and
% value functions for SDF


resfile = 'res_20191112_bench.mat';

load(resfile);

TOL=1e-4;
MAXIT=2500;

%% Initialize
betaS_g = mobj.Params.betaS;
betaI_g = mobj.Params.betaI;
betaB_g = mobj.Params.betaB;

% Intial guess for LT bond price
qT_ss = betaS_g / (1 - betaS_g * mobj.Params.delta);

% Initial guess for cB stream price
qCB_ss = betaB_g / (1 - betaB_g) * stv{1}.Sol.cB;
% Initial guess for cI stream price
qCI_ss = betaI_g / (1 - betaI_g) * stv{1}.Sol.cI;
% Initial guess for cS stream price
qCS_ss = betaB_g / (1 - betaS_g) * stv{1}.Sol.cS;
pricenames = {'qT','qCB','qCI','qCS'};


npt = size(mobj.Vfct.Vals,1);

Pricefct = grid.LinearInterpFunction( mobj.Vfct.SSGrid, ...
	repmat( [qT_ss, qCB_ss, qCI_ss, qCS_ss], npt, 1 ) );

%% Compute SDFs (once)
SDFS = zeros(npt,mobj.Exogenv.exnpt);
SDFB = zeros(npt,mobj.Exogenv.exnpt);
SDFI = zeros(npt,mobj.Exogenv.exnpt);
transtens = zeros(1+mobj.NSTEN,mobj.Exogenv.exnpt,npt);

pointmat = mobj.Vfct.SSGrid.Pointmat;
polmat = mobj.Pfct.Vals;
transmat=mobj.Tfct.Vals;
EXNPT=mobj.Exogenv.exnpt;
NSTEN=mobj.NSTEN;

parfor i=1:npt
	pt=pointmat(i,:);
	trans=transmat(i,:);%trans=mobj.Tfct.Vals(i,:);
	solvec=polmat(i,:);
	cB=exp(solvec(6));
	cI=exp(solvec(7));
	cS=exp(solvec(8));
	p0=exp(solvec(4));
	transpts=reshape([(1:EXNPT),trans],EXNPT,NSTEN+1);%transpts=reshape([(1:mobj.Exogenv.exnpt),trans],mobj.Exogenv.exnpt,mobj.NSTEN+1);
	vtrans=mobj.evaluateVal(transpts)';%vtrans=mobj.evaluateVal(transpts)';
	expStruct=mobj.computeExpectations(pt(1),transpts,p0,vtrans);
   
	% Pull SDF of a particular agent (here, saver)
    cvec=[cB,cI,cS];
    hvec=[mobj.Params.KBbar,mobj.Params.HI,mobj.Params.HS];
    xi=mobj.Exogenv.pts_perm(pt(1),3);
    U_vec=cvec.^(1-xi) .* hvec.^xi;
    U_norm=U_vec.^(1-1/mobj.Params.psi)./cvec;
	SDFB(i,:)=expStruct.Conditional.SDFB(:,1)' / U_norm(1);
	SDFI(i,:)=expStruct.Conditional.SDFI(:,1)' / U_norm(2);
	SDFS(i,:)=expStruct.Conditional.SDFS(:,1)' / U_norm(3);
   
	% Store transition states for interpolating future payoffs in the
	% iteration part later
	transtens(:,:,i) = transpts';
end

%% Iterate on payoffs
iter=1;
mtrans = mobj.Exogenv.mtrans;

params = mobj.Params;

while true
	resmat = Pricefct.Vals;
	parfor i=1:npt
		
        qvec=zeros(1,4);
		prnext = mtrans(pointmat(i,1),:);
        SDFmat = [SDFB(i,:)',SDFI(i,:)',SDFS(i,:)'];
		futurePrice = Pricefct.evaluateAt( transtens(:,:,i)' )';
		
		% LT Bond
		payoff = 1 + params.delta * futurePrice(:,1);
		qT = prnext * ( SDFmat(:,3) .* payoff );
		qvec(1) = qT;
        
		% Consumptions
		pol = polmat(i,:);
		
        for a=0:2
    		c = exp( pol(6+a) );
        	payoff = c + futurePrice(:,2+a);
            qC = prnext * ( SDFmat(:,1+a) .* payoff );
            qvec(2+a) = qC;
        end
        
		resmat(i,:)=qvec;
	end
	
	[dist,wh]=max(abs(resmat(:)-Pricefct.Vals(:)));     
	[wh_1,wh_2]=ind2sub(size(resmat),wh);
	disp(['-- Iteration: ',num2str(iter),', max distance: ',num2str(dist),' in ',char(pricenames(wh_2)), ...
		' at point ',num2str(wh_1),': ',num2str(mobj.Vfct.SSGrid.Pointmat(wh_1,:))]);

	Pricefct=Pricefct.fitTo(resmat);
	
	if dist <= TOL
		disp('Converged!');
		break;
	elseif iter>=MAXIT
        disp('Max.iter. exceeded.');
        break;
	end
	
	iter=iter+1;
end

%% Apply to simulation
load(['sim_',resfile]);
nvar=size(simseries,2);
tmp_statevec = statevec(2:end);
tmp_simseries = simseries(:,4:8);
tmp = Pricefct.evaluateAt( [ tmp_statevec , tmp_simseries ] )';
qT = tmp(:,1);
qCB = tmp(:,2);
qCI = tmp(:,3);
qCS = tmp(:,4);

rT = log(1./qT+mobj.Params.delta);

simseries=[simseries,qT,rT,qCB,qCI,qCS];
varnames=[varnames,{'qT','rT','qCB','qCI','qCS'}];

mobj=[];
save(['pzns_',resfile]);


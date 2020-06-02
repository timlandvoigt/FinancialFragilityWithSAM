classdef SAMIntermediaryModel < model.DSGEModel
    
    properties (SetAccess=protected)
        NSTEN % number of endogenous state variables: Vfct.Ndim-1 
        NSTEX % number of exogenous state variables
        NSOL % number of solution vars
        NV % number of forecasting variables
        NADD % number of additional endogenous variables
        NCOND % number of conditional expectations of period-ahead variables        
        Sol_names % NSOLx1 cell array with names of solution variables
                   % must be in order of functions of Pfct 
        V_names % NVx1 cell array with names of forecasting variables
                   % must be in order of functions of Vfct      
        Sol_baseguess           
        En_names % NSTENx1 cell array with names of endog state variables           
        Ex_names % NSTEXx1 cell array with names of exog state variables           
        Add_names % NADDx1 cell array with names of additional vars
        Cond_names %NCONDx1 cell array with names of conditional variables        
        Params % structure array with parameters
        Exogenv % structure array with fields mtrans, pts_perm, pts_all
                % specifying exog. Markov processes
        Vfct % ApproxFunction object for iteration-relevant functions
        Pfct % ApproxFunction object for solution jump variables
        Tfct % ApproxFunction object for transition of state variable(s) (optional)
        Tempfct % ApproxFunction for intermediate computations 
        Basegrid % base grid 
        Jacobian % Function evaluating the nonlinear system Jacobian
    end    
    
    
    methods
        % constructor
        function obj=SAMIntermediaryModel(params,endogenv,exogenv,vfct,pfct,tfct,tempfct,basegrid,jac)
            % call superclass constructor
            obj=obj@model.DSGEModel(params,endogenv,exogenv,vfct,pfct,tfct);
            obj.Tempfct=tempfct;
            obj.Basegrid=basegrid; 
            obj.NCOND=length(endogenv.condnames);          
            obj.Cond_names=reshape(endogenv.condnames,1,obj.NCOND);
            obj.Jacobian=jac;
        end
        
        function [nextst,outstr]=calcStateTransition(obj,point,solvec,mode,varargin)
            params=obj.Params;
            % unpack relevant params
            shareB=params.shareB;
            shareI=params.shareI;
            shareS=params.shareS;
            tau=params.tau;
            delta=params.delta;
            phiK=params.phiK;
            nuREO=params.nuREO;
            KBbar=params.KBbar;
            %xi=params.xi;
            SREO=params.SREO;
            iota_om=params.iota_om;
            HI=params.HI;
            HS=params.HS;
            sigma_eps=params.sigma_eps;
            zeta=params.zeta;
            xiL=params.xiL;
            pi=params.pi;
            pgr_max=params.pgr_max;
            zetaombar=params.pgr_max;
            pgr_min=params.pgr_min;
            mu_G=params.mu_G;
            iota_p=params.iota_p;
            indexMB=params.indexMB;
            indexAB=params.indexAB;
            indexTail=params.indexTail;
            tauL=params.tauL;
            theta=params.theta;
            liq_thresh=params.liq_thresh;
            etaB=params.etaB;
            if isfield(params,'rho_om')
                rho_om=params.rho_om;
            else
                rho_om=0.95;
            end
            
            
            % extract state variables
            exst=point(1);
            AB=point(2);
            MB=point(3);
            KREO=point(4);
            WI=point(5);
            BG=point(6);
            Y=obj.Exogenv.pts_perm(exst,1);
            sig2_om=obj.Exogenv.pts_perm(exst,2);
            xi=obj.Exogenv.pts_perm(exst,3);
            KB=KBbar-KREO;
            
            fullendogvars=[AB,MB,KREO,WI,BG];            
            
            % compute next period's states
            if isempty(varargin)
                State_next=obj.evaluateTrans(point);
                thisval=obj.evaluateVal(point);
            else
                State_next=varargin{1};
                thisval=varargin{2};
            end
            
            OmA=thisval(7);
            OmM=thisval(8);
            CscrI=thisval(15);
            
            % extract solution variables
            qA=exp(solvec(1));
            qM=exp(solvec(2));
            qf=exp(solvec(3));
            p=exp(solvec(4));
            pREO=exp(solvec(5));
            cB=exp(solvec(6));
            cI=exp(solvec(7));
            cS=exp(solvec(8));
            ZR=exp(solvec(9));
            MS_g=solvec(10);
            lamLTV=solvec(11);
            lamLTVplus=max(0,lamLTV)^3;
            
            % intraperiod FOC for current interest rate
            rstar=(1-qM)/qA;
            
            % Maintenance technology
            if params.maint_by_price
                nuK = params.nuK;
                nuREO = params.nuREO;
            else
                nuK = params.nuK_p / p;
                nuREO = params.nuREO_p / p;
            end
            if params.rebate_maint
                nuK_budget = 0;
            else
                nuK_budget = nuK;
            end                
                
            
            % Quadrature way
            [nodes, weights] = GaussHermite(11);
            
            QK = (1-nuK - (1-ZR)*phiK*lamLTVplus)*p*KB;
            QA = ( (1-tau)*AB + delta*(1-ZR)*(OmA*AB) );
            QM = ( (1-delta + delta*ZR)*MB + delta*(1-ZR)*(OmM*MB) );
            
%             [ZN_old, ZK_old, ZM_old, ZA_old, jacout_old] = SAMIntermediaryModel.Mpayoff_logn_quad(QK, QA, QM, ...
%                 sig2_om, rho_om, iota_om, zetaombar, indexAB, indexMB, ...
%                 nodes, weights, true);
            
            % LIQUIDITY
                
            % Input values for testing
            QK_liq = (1 - nuK) * liq_thresh * p * KB;
            
            [ZN, ZK, ZM, ZA, ...
                ~, ~, ~, ~, ...
                ~, ~, ~, ~, ...
                ~, ~, ~, ~, ...
                Z_vals, jacout] = SAMIntermediaryModel.Mpayoff_logn_quad_liq(QK, QA, QM, QK_liq, ...
                KB, MB, AB, p, sig2_om, rho_om, iota_om, zetaombar, ...
                indexAB, indexMB, indexTail,theta, etaB,...
                nodes, weights, true);
            
            % END CHANGE
            
            
            % intermediary bankruptcy and bailout
            VI=WI+CscrI;
            if params.noBR
                F_eps=1;
                F_eps_cond=0;
            else
                F_eps=normcdf(VI,0,sigma_eps);
                F_eps_cond=normpdf(VI/sigma_eps);
            end
            F_eps_minus=-sigma_eps*F_eps_cond;
            F_eps_plus=sigma_eps*F_eps_cond;
           
            % government debt and spending
            bailout=F_eps_plus - (1-F_eps)*(WI - zeta*delta*(1-ZR)*(ZA*qA*AB+ZM*qM*MB));
            Tvec=[shareB, shareI, shareS]' * (tauL*(BG + bailout));
            BGnext=(1-tauL)*(BG + bailout)/qf; % gov debt EOP

            G=tau*Y - tau*(ZA*AB);
            
            
            % borrower budget
            rent=xi*cB/((1-xi)*KBbar);            
            Mstar=(cB - ( (1-tau)*shareB*Y - ZA*(1-tau)*AB - ZM*MB*((1-delta) + delta*ZR) ...
                - p*SREO*KREO - rent*KREO - nuK_budget.*p*ZK*KB - Tvec(1)) )/(ZR*ZN);
            AB_g=ZR*ZN*rstar*Mstar+delta*(1-ZR)*ZA*AB;
            MB_g=ZR*ZN*Mstar+delta*(1-ZR)*ZM*MB;
            Kstar=(SREO*KREO + ZR*ZK*KB)/(ZR*ZN);
            
            % mortgage market clearing
            AS_g = AB_g*MS_g/MB_g;
            MI_g = MB_g - MS_g;
            AI_g = AB_g - AS_g;
            
            % intermediary and saver budget
            X = (1-ZK)*KB*(pREO-nuREO*p)/MB;
            debtB = (X + ZM*(1-delta + delta*ZR))*MB + ZA*AB + delta*(1-ZR)*(ZM*qM*MB+ZA*qA*AB);
            wealthS = debtB - WI + BG; % total debt held by savers BOP
            BSnext = -(cS - ((1-tau)*shareS*Y + wealthS - qM*MS_g - qA*AS_g - nuK_budget.*p*HS - Tvec(3)))/qf;  % saver debt EOP
            cy = xiL*cS/((1-xi-xiL)*BSnext);
            if cy<0
                cy=1e20;
            end
            % impose MC in mortgage market to compute intermediary budget
            KREO_g = (1-SREO)*KREO + (1-ZK)*KB;
            BInext=(cI - ( F_eps*WI - F_eps_minus + (1-tau)*shareI*Y - qM*MI_g - qA*AI_g ...
                + (SREO-nuREO)*p*KREO + rent*KREO - pREO*(1-ZK)*KB - nuK_budget.*p*HI - Tvec(2)))/qf;
                                    
            %BInext = BSnext - BGnext;  % deposits EOP (market clearing)
            
            
            if mode>0
                % simulation, mode contains number of next period's state
                exst=mode;
                exnpt=size(obj.Exogenv.pts_all,1);
                cind=obj.Exogenv.pts_all(exst,end);
                sig2_om_next=obj.Exogenv.pts_all(exst,2);
                ABnext=State_next(exst,1);
                MBnext=State_next(exnpt+exst,1);
                KREOnext=State_next(2*exnpt+exst,1);
                WInext=State_next(3*exnpt+exst,1);
                BGtrans=State_next(4*exnpt+exst,1);
                nextst=[cind,ABnext,MBnext,KREOnext,WInext,BGtrans];
                expstr=obj.computeExpectations(exst, nextst, p);
                p_next = expstr.p_next;
                Mpayoff_next = expstr.Mpayoff_next;
                Apayoff_next = expstr.Apayoff_next;
                pgr=p_next/p;
                pgr(pgr>pgr_max)=pgr_max;
                pgr(pgr<pgr_min)=pgr_min;
                if indexTail == 1
                    pgr = pgr + (1-pgr_max);
                end
                zeta_next = (1 - iota_p + iota_p * pgr) ./ mu_G;
                ABtrans = AB_g * (indexAB*zeta_next + 1-indexAB) / pi;
                MBtrans = MB_g * (indexMB*zeta_next + 1-indexMB) / pi;
                KREOtrans = KREO_g ./ mu_G;
                WItrans = (Mpayoff_next.*MBtrans + Apayoff_next.*ABtrans)*(1-MS_g/MB_g) - BInext/(mu_G*pi);
                BGtrans = BGnext/(mu_G*pi);
                nextst=[cind,ABtrans,MBtrans,KREOtrans,WItrans,BGtrans];
                
            else
                % solution mode, compute next period's state variables for all
                % possible Markov states
                sig2_om_next=obj.Exogenv.pts_all(:,2);
                cind=obj.Exogenv.pts_all(:,end);
                nst=obj.Exogenv.exnpt;
                ABnext=State_next(1:nst,1);
                MBnext=State_next(nst+(1:nst),1);
                KREOnext=State_next(2*nst+(1:nst),1);
                WInext=State_next(3*nst+(1:nst),1);
                BGtrans=State_next(4*nst+(1:nst),1);
                % matrix of next period states
                nextst=[cind,ABnext,MBnext,KREOnext,WInext,BGtrans];
            end
            
            nn = 1;
            
            ZN_str = Z_vals(nn); nn = nn + 1;
            ZK_str = Z_vals(nn); nn = nn + 1;
            ZM_str = Z_vals(nn); nn = nn + 1;
            ZA_str = Z_vals(nn); nn = nn + 1;
            
            ZN_liq = Z_vals(nn); nn = nn + 1;
            ZK_liq = Z_vals(nn); nn = nn + 1;
            ZM_liq = Z_vals(nn); nn = nn + 1;
            ZA_liq = Z_vals(nn); nn = nn + 1;
            
            addvars=struct(...
                'ZK',ZK,...
                'ZN',ZN,...
                'ZA',ZA,...
                'ZM',ZM,...
        		'ZK_str',ZK_str,...
                'ZN_str',ZN_str,...
                'ZA_str',ZA_str,...
                'ZM_str',ZM_str,...
                'ZK_liq',ZK_liq,...
                'ZN_liq',ZN_liq,...
                'ZA_liq',ZA_liq,...
                'ZM_liq',ZM_liq,...
                'F_eps',F_eps,...
                'F_eps_minus',F_eps_minus,...
                'F_eps_plus',F_eps_plus,...
                'VI', VI, ...
                'Mstar',Mstar,...
                'Lstar',ZR*ZN*Mstar,...
                'Kstar',Kstar,...
                'AB_g',AB_g,...
                'MB_g',MB_g,...
                'BSnext',BSnext,...
                'BInext',BInext,...
                'BGnext',BGnext,...
                'cy',cy,...
                'bailout',bailout,...
                'X',X,...
                'KREO_g',KREO_g,...
                'KB',KB,...
                'MB',MB,...
                'AB',AB,...
                'OmA',OmA,...
                'OmM',OmM,...
                'sig2_om', sig2_om,...
                'sig2_om_next',sig2_om_next,...
                'rent',rent,...
                'G',G);
            
            outstr=struct;
            outstr.addvars=addvars;
            outstr.exstvec=[Y;sig2_om;xi];
            outstr.Jvec=[F_eps,F_eps_cond,jacout];
            outstr.fullendogvars=fullendogvars;
            
        end
        
        
        function [fx,J,V]=calcEquations(obj,exst,nextst,solvec,instr,mode,varargin)
            % allocate result
            fx=zeros(obj.NSOL,1);
            J=[];
            
            % unpack params
            params=obj.Params;
            betaB=params.betaB;
            betaI=params.betaI;
            betaS=params.betaS;
            psi=params.psi;
            %xi=params.xi;
            mu_G=params.mu_G;
            phiK=params.phiK;
            phiI=params.phiI;
            tau=params.tau;
            xiL=params.xiL;
            delta=params.delta;
            iota_om=params.iota_om;
            iota_p=params.iota_p;
            HI=params.HI;
            HS=params.HS;
            KBbar=params.KBbar;
            pi=params.pi;
            mu_kappa=params.mu_kappa;
            sig_kappa=params.sig_kappa;
            sigma_eps=params.sigma_eps;
            pgr_max=params.pgr_max;
            pgr_min=params.pgr_min;
            indexMB=params.indexMB;
            indexAB=params.indexAB;
            indexTail=params.indexTail;
            internalize_liq=params.internalize_liq;
            psiS0=params.psiS0;
            psiS1=params.psiS1;
            
            % extract endogeous variables 
            qA=exp(solvec(1));
            qM=exp(solvec(2));            
            qf=exp(solvec(3));            
            p=exp(solvec(4));
            pREO=exp(solvec(5));
            cB=exp(solvec(6));
            cI=exp(solvec(7));  
            cS=exp(solvec(8));  
            ZR=exp(solvec(9));   
            MS_g=solvec(10);
            lamLTV=solvec(11);
            lamI=solvec(12);
            muS=solvec(13);
            lamS=solvec(14);

            % multiplier transformations
            lamLTVplus=max(0,lamLTV)^3;
            lamLTVminus=max(0,-lamLTV)^3;           
            lamIplus=max(0,lamI)^3;
            lamIminus=max(0,-lamI)^3;           
            muSplus=max(0,muS)^3;
            muSminus=max(0,-muS)^3;           
            lamSplus=max(0,lamS)^3;
            lamSminus=max(0,-lamS)^3;           
            % intraperiod FOC for current interest rate
            rstar=(1-qM)/qA;
            
            % extract some other state-dependent variables
            envec=instr.addvars;
            ZK=envec.ZK;
            ZA=envec.ZA;
            ZM=envec.ZM;
            ZN=envec.ZN;
            ZK_str=envec.ZK_str;
            ZA_str=envec.ZA_str;
            ZM_str=envec.ZM_str;
            ZN_str=envec.ZN_str;
            ZK_liq=envec.ZK_liq;
            ZA_liq=envec.ZA_liq;
            ZM_liq=envec.ZM_liq;
            ZN_liq=envec.ZN_liq;
            AB_g=envec.AB_g;    
            MB_g=envec.MB_g;
            BSnext=envec.BSnext;
            BInext=envec.BInext;
            BGnext=envec.BGnext;
            cy=envec.cy;
            KB=envec.KB;
            AB=envec.AB;
            MB=envec.MB;
            OmA=envec.OmA;
            OmM=envec.OmM;
            Mstar=envec.Mstar;
            Kstar=envec.Kstar;
            KREO_g=envec.KREO_g;
            VI=envec.VI;
            xi=instr.exstvec(3);
            
            % Maintenance technology
            if params.maint_by_price
                nuK = params.nuK;
            else
                nuK = params.nuK_p / p;
            end
                                  
            % probabilities and states to compute expectation terms
            prnext=obj.Exogenv.mtrans(exst,:);          
            
            % projection evaluation
            if nargin==7
                % Exactly 7 arguments were passed.
                % Means expectations have already been computed (fast
                % solution mode)
                expStruct = varargin{1};
            elseif nargin==6
                % Nothing passed in ((conditional moments and EE errors in
                % simulation)
                expStruct = obj.computeExpectations(exst,nextst,p);
            else
                % trans, vtrans, and vtrans_def were passed in (legacy)  
                expStruct = obj.computeExpectations(exst,nextst,p,varargin{:});
            end            

            exp_OmA_B = expStruct.exp_OmA_B;
            exp_OmM_B = expStruct.exp_OmM_B;
            exp_OmK_B = expStruct.exp_OmK_B;
            exp_OmA_I = expStruct.exp_OmA_I;
            exp_OmM_I = expStruct.exp_OmM_I;
            exp_OmMA_S = expStruct.exp_OmMA_S;
            exp_SDFI = expStruct.exp_SDFI;
            exp_SDFS = expStruct.exp_SDFS;
            exp_pREO_next = expStruct.exp_pREO_next;
			exp_CscrI_next = expStruct.exp_CscrI_next;
            p_next = expStruct.p_next;
            Mpayoff_next = expStruct.Mpayoff_next;
            Apayoff_next = expStruct.Apayoff_next;
            CE_vec = expStruct.CE_vec;
            
             
            % compute certainty equiv. and SDFs
            V_vec=zeros(3,1);
            U_norm=zeros(3,1); % to normalize SDFs back
            cvec=[cB,cI,cS];
            hvec=[KBbar,HI,HS];
            beta_vec=[betaB,betaI,betaS];
            xiLvec=[0,0,xiL];
            U_vec=cvec.^(1-xi-xiLvec) .* hvec.^xi;
            U_vec(3)=U_vec(3)*BSnext^xiL;
            psiinv_vec=1/psi*ones(1,3);
            for i=1:3
                if psiinv_vec(i)==1
                   V_vec(i)=U_vec(i)^(1-beta_vec(i)) *  CE_vec(i)^beta_vec(i); 
                else
                   V_vec(i)=( (1-beta_vec(i))*U_vec(i)^(1-psiinv_vec(i)) + beta_vec(i)* CE_vec(i)^(1-psiinv_vec(i)) )^(1/(1-psiinv_vec(i))); 
                end
                U_norm(i)=U_vec(i)^(1-psiinv_vec(i))/cvec(i);  
            end      
            
            if mode==3
                Vnext=zeros(obj.Vfct.Nof,1);
                Vnext(4:6)=V_vec;
                V{1}=Vnext;
                return;
            end
                       
            U_norm_p=U_norm;
            exp_OmM_B_norm=exp_OmM_B/U_norm_p(1);
            exp_OmA_B_norm=exp_OmA_B/U_norm_p(1);
            exp_OmK_B_norm=exp_OmK_B/U_norm(1);  
            rbar=AB/MB;
            kappabar = (1-exp_OmM_B_norm-rbar*exp_OmA_B_norm)*(1-delta*ZM*MB/(ZN*Mstar)) ...
                             - exp_OmA_B_norm*(rstar-rbar) - p*phiK*lamLTVplus*(ZN*Kstar-ZK*KB)/(ZN*Mstar);
            % borrower FOCs
            fx(1) = 1-lamLTVplus - exp_OmM_B_norm - rstar*exp_OmA_B_norm;
            fx(2) = p*(1-lamLTVplus*phiK) - exp_OmK_B_norm;
            fx(3) = ZR - obj.Refi(kappabar,mu_kappa,sig_kappa);
            
            exp_OmA_I_norm=exp_OmA_I/U_norm_p(2);
            exp_OmM_I_norm=exp_OmM_I/U_norm_p(2);
            exp_pREO_next_norm=exp_pREO_next/U_norm(2);
            % lender FOCs
            fx(4)=qA*(1-phiI*lamIplus) - exp_OmA_I_norm;
            fx(5)=qM*(1-phiI*lamIplus) - exp_OmM_I_norm;
            fx(6)=pREO - exp_pREO_next_norm;
            fx(7)=qf - lamIplus - exp_SDFI/U_norm(2);
            
            % saver FOC
            rbar_g=AB_g/MB_g;
            fx(8)=qf - cy - exp_SDFS/U_norm(3) - lamSplus;
            fx(9)=qM + rbar_g*qA + psiS0*(MS_g)^(psiS1-1) - exp_OmMA_S/U_norm(3) - muSplus;
            
            % borrower constraint
            fx(10)=phiK*p*Kstar - Mstar - lamLTVminus; 

            % intermediary constraint
            fx(11)=phiI*(1-MS_g/MB_g)*(qA*AB_g+qM*MB_g) - BInext - lamIminus;   
            
            % saver no-shorting constraint
            fx(12)=MS_g - muSminus;
            fx(13)=BSnext - lamSminus;
            
            % bond market clearing
            fx(14)=BInext + BGnext - BSnext;
                        
            % Transitions
            nst = obj.Exogenv.exnpt;
            pgr=p_next/p;
            pgr(pgr>pgr_max)=pgr_max;
            pgr(pgr<pgr_min)=pgr_min;
            if indexTail == 1
                pgr = pgr + (1-pgr_max);
            end
            zeta_next = (1 - iota_p + iota_p * pgr) ./ mu_G;
            ABtrans = AB_g * (indexAB*zeta_next + 1-indexAB) / pi;
            MBtrans = MB_g * (indexMB*zeta_next + 1-indexMB) / pi;          
            KREOtrans = ones(nst,1) * KREO_g ./ mu_G;
            WItrans = (Mpayoff_next.*MBtrans + Apayoff_next.*ABtrans)*(1-MS_g/MB_g) - BInext/(mu_G*pi);
            BGtrans = ones(nst,1) * BGnext/(mu_G*pi);
            
            % Jacobian computation (numeric)
            if nargout > 1 && mode == 0 && ~isempty(obj.Jacobian)
                fullendogvars=instr.fullendogvars;
                Y=instr.exstvec(1);
                sig2_om=instr.exstvec(2);
                
                % Compute additional things needed for Jacobian:
                F_eps = instr.Jvec(1);
                F_eps_cond = instr.Jvec(2);
                ZZN_str = ZN_str;
                ZZK_str = ZK_str;
                ZZA_str = ZA_str;
                ZZM_str = ZM_str;
                ZZN_liq = ZN_liq;
                ZZK_liq = ZK_liq;
                ZZA_liq = ZA_liq;
                ZZM_liq = ZM_liq;
                dZZN_str = instr.Jvec(3);
                dZZK_str = instr.Jvec(4);                
                dZZA_str = instr.Jvec(5);
                dZZM_str = instr.Jvec(6);   
                dZZN_liq = instr.Jvec(7);
                dZZK_liq = instr.Jvec(8);                
                dZZA_liq = instr.Jvec(9);
                dZZM_liq = instr.Jvec(10); 
                dF_eps = normpdf(VI,0,sigma_eps);
                dF_eps_cond = dF_eps*VI;
                nodes = GaussHermite(11);
                
                % Compute Jacobian
                                
                % Needs to have the same order as constructJacobian
                args = [{Y,sig2_om,xi},num2cell(fullendogvars),{qA,qM,qf,p,pREO,cB,cI,cS,ZR,MS_g,...
                lamLTVplus,lamLTVminus,lamIplus,lamIminus,muSplus,muSminus,lamSplus,lamSminus,...
                ZZN_str,ZZK_str,ZZA_str,ZZM_str,dZZN_str,dZZK_str,dZZA_str,dZZM_str,...
                ZZN_liq,ZZK_liq,ZZA_liq,ZZM_liq,dZZN_liq,dZZK_liq,dZZA_liq,dZZM_liq,...
                F_eps, F_eps_cond, dF_eps, dF_eps_cond,...				
                exp_OmA_B,exp_OmM_B,exp_OmK_B,exp_OmA_I,exp_OmM_I,exp_OmMA_S,exp_SDFI,exp_SDFS,exp_pREO_next,...
                OmA, OmM, nodes(6)}];

                % Jacobian with respect to actual solution variables
                tmp_Jacobian = obj.Jacobian(args{:});

                nsol = size(solvec,1);
                J = zeros(nsol);

                % Now use chain rule to compute Jacobian with respect to
                % solvec
                
                % Exponentiated variables
                J(:,1:9) = tmp_Jacobian(1:end,1:9) .* repmat( exp(solvec(1:9))' ,nsol,1);          
                J(:,10) = tmp_Jacobian(1:end,10);

                % Multipliers
                start_idx = 10;
                
                for mult_idx = 1:nsol-start_idx
                    mult = solvec(start_idx + mult_idx);
                    if mult > 0
                        J(:,start_idx + mult_idx) = ...
                            tmp_Jacobian(1:end,start_idx + 2*(mult_idx-1) + 1) .* repmat( 3*mult^2 ,nsol,1);
                    else
                        J(:,start_idx + mult_idx) = ...
                            tmp_Jacobian(1:end,start_idx + 2*(mult_idx-1) + 2) .* repmat( -3*mult^2 ,nsol,1);
                    end
                end
            end
            
            V=cell(3,1);
            % marginal value functions
            if mode==1
                % Output new values for time iteration
                % risk-taker value function evaluated at zero risk-taker
                Vnext=zeros(obj.Vfct.Nof,1);
                Vnext(1)=cB;
                Vnext(2)=cI;
                Vnext(3)=cS;
                Vnext(4)=V_vec(1);
                Vnext(5)=V_vec(2);
                Vnext(6)=V_vec(3);
                Vnext(7)=exp_OmA_B_norm;
                Vnext(8)=exp_OmM_B_norm;
                Vnext(9)=exp_OmK_B_norm;
                Vnext(10)=qA;
                Vnext(11)=qM;
                Vnext(12)=ZR;
                Vnext(13)=p;
                Vnext(14)=pREO;
				Vnext(15)=exp_CscrI_next/U_norm(2);
                Vnext(16)=lamLTV;
                Vnext(17)=BSnext;
                V{1}=Vnext;
                % state transition
                V{2}=[ABtrans; MBtrans; KREOtrans; WItrans; BGtrans]';
                if ~isreal(Vnext)
                    disp(dang);
                end
                
        
            elseif mode==2
                
                % Evaluation during simulation. Output conditional
                % variables
                
                % SDFs
                SDFB = expStruct.Conditional.SDFB/U_norm(1);
                SDFI = expStruct.Conditional.SDFI/U_norm(2);
                Apayoff_next = expStruct.Conditional.Apayoff_next;
                Mpayoff_next = expStruct.Conditional.Mpayoff_next;
                Del_BC_K_next = expStruct.Conditional.Del_BC_K_next;
                Del_BC_M_next = expStruct.Conditional.Del_BC_M_next;
               
                SDF.SDFI = SDFI;
                SDF.SDFB = SDFB;    
                
                % Define returns to intermediaries on IO and PO strips
                Apayoff_next_index=(indexAB.*zeta_next + 1-indexAB)/pi .* Apayoff_next;
                Mpayoff_next_index=(indexMB.*zeta_next + 1-indexMB)/pi .* Mpayoff_next;
                retA = Apayoff_next_index/qA;
                retM = Mpayoff_next_index/qM;
                retT = rstar*Apayoff_next_index + Mpayoff_next_index;
                expRP_A = prnext * retA;
                expRP_M = prnext * retM;
                expRP_T = prnext * retT;
                CovP_A = prnext* ( (SDF.SDFI - prnext*SDF.SDFI) .* (retA - expRP_A) );


                % New version (DG)
                Del_BC_K = internalize_liq * prnext * (SDFB.*Del_BC_K_next);
                Del_BC_M = internalize_liq * prnext * (SDFB.*Del_BC_M_next);
                
                condvars = struct('expRP_A',expRP_A, ...
                                  'expRP_M',expRP_M,...
                                  'expRP_T',expRP_T,...
                                  'CovP_A',CovP_A,...
                                  'Del_BC_K',Del_BC_K,...
                                  'Del_BC_M',Del_BC_M);                        
                
                Wtrans.AB = ABtrans;                          
                Wtrans.MB = MBtrans;
                Wtrans.KREO = KREOtrans;
                Wtrans.WI = WItrans;
                Wtrans.BG = BGtrans;
                V = {condvars,Wtrans,SDF};
            end
        end
    
        
        function expStruct = computeExpectations( obj, exst, nextst, p0, varargin )
            % This function computes E[ ] terms in equilibrium conditions
            
            % unpack params
            params=obj.Params;
            betaB=params.betaB;
            betaI=params.betaI;
            betaS=params.betaS;
            sigmaB=params.sigmaB;
            sigmaI=params.sigmaI;
            sigmaS=params.sigmaS;
            psi=params.psi;
            pi=params.pi;
            %xi=params.xi;
            mu_G=params.mu_G;
            nuREO=params.nuREO;
            SREO=params.SREO;
            tau=params.tau;
            delta=params.delta;
            phiK=params.phiK;
            iota_om=params.iota_om;
            iota_p=params.iota_p;
            HI=params.HI;
            HS=params.HS;
            KBbar=params.KBbar;
            xiL=params.xiL;
            pi=params.pi;
			sigma_eps=params.sigma_eps;
            pgr_max=params.pgr_max;
            pgr_min=params.pgr_min;
            zetaombar=pgr_max;
            indexMB=params.indexMB;
            indexAB=params.indexAB;
            indexTail=params.indexTail;
            theta=params.theta;
            liq_thresh=params.liq_thresh;
            etaB=params.etaB;
            internalize_liq=params.internalize_liq;
            shareB=params.shareB;
            tauL=params.tauL;
            if isfield(params,'rho_om')
                rho_om=params.rho_om;
            else
                rho_om=0.95;
            end
            
            % probabilities and states to compute expectation terms
            prnext=obj.Exogenv.mtrans(exst,:);   
            Y_next=obj.Exogenv.pts_perm(nextst(:,1),1);
            sig2_om_next=obj.Exogenv.pts_perm(nextst(:,1),2);
            xi_next=obj.Exogenv.pts_perm(nextst(:,1),3);
            
            % projection evaluation
            onlyCE=false;
            if nargin>4
                Pol_next=varargin{1};
                if nargin>5
                    expStruct=varargin{2};
                    onlyCE=true;
                end
            else
                Pol_next=obj.evaluateVal(nextst)';
                if size(Pol_next,1)==1
                    prnext=1;
                end
            end
          
            CB_next=Pol_next(:,1);
            CI_next=Pol_next(:,2);
            CS_next=Pol_next(:,3);
            VB_next=Pol_next(:,4);
            VIHH_next=Pol_next(:,5);
            VS_next=Pol_next(:,6);
            exp_OmA_B_next=Pol_next(:,7);
            exp_OmM_B_next=Pol_next(:,8);
            exp_OmK_B_next=Pol_next(:,9);
            qA_next=Pol_next(:,10);
            qM_next=Pol_next(:,11);
            ZR_next=Pol_next(:,12);
            p_next=Pol_next(:,13);
            pREO_next=Pol_next(:,14);
			CscrI_next=Pol_next(:,15);
            lamLTV_next=Pol_next(:,16);
            lamLTVplus_next=max([lamLTV_next,zeros(size(lamLTV_next))],[],2).^3;
			BSnext_next=Pol_next(:,17);

            
            % Maintenance technology
            if params.maint_by_price
                nuK = params.nuK;
                nuREO = params.nuREO;
            else
                nuK = params.nuK_p ./ p_next;
                nuREO = params.nuREO_p ./ p_next;
            end
                        
            % intermediary bankruptcy and bailout
            VI_next=nextst(:,5)+CscrI_next;
            if params.noBR
                F_eps_next=1;
                F_eps_cond_next=0;
            else
                F_eps_next=normcdf(VI_next,0,sigma_eps);
                F_eps_cond_next=normpdf(VI_next./sigma_eps);
            end
            F_eps_minus_next=-sigma_eps.*F_eps_cond_next;
            F_eps_plus_next=sigma_eps.*F_eps_cond_next;
            
            % compute certainty equiv. and SDFs for borrower and saver
            UB_next=CB_next.^(1-xi_next).*KBbar.^xi_next;
            UI_next=CI_next.^(1-xi_next).*HI.^xi_next;
            US_next=CS_next.^(1-xi_next-xiL) .* HS.^xi_next .* BSnext_next.^xiL;
            Vnext_mat={VB_next,VIHH_next,VS_next};
            Cnext_mat={CB_next,CI_next,CS_next};
            Unext_mat={UB_next,UI_next,US_next};
            sigma_vec=[sigmaB,sigmaI,sigmaS];
            beta_vec=[betaB,betaI,betaS];
            psiinv_vec=1/psi*ones(1,3);
            SDF_mat=cell(3,1);
            for i=1:3
                % certainty equivalent
                Vnext_tmp=Vnext_mat{i};
                Cnext_tmp=Cnext_mat{i};
                Unext_tmp=Unext_mat{i};
                if sigma_vec(i)==1          
                    % log
                    CE_tmp=exp(prnext * log(mu_G.*Vnext_tmp));
                else                
                    CE_tmp=(prnext * (mu_G.*Vnext_tmp).^(1-sigma_vec(i)))^(1/(1-sigma_vec(i)));
                end
                Unext_tmp=Unext_tmp.^(1-psiinv_vec(i))./Cnext_tmp;
                SDF_tmp=beta_vec(i)* mu_G.^(-sigma_vec(i))...
                    .* Unext_tmp .* (Vnext_tmp/CE_tmp).^(psiinv_vec(i)-sigma_vec(i));
                CE_vec(i)=CE_tmp;
                SDF_mat{i}=SDF_tmp;
            end         
            
            if onlyCE
                expStruct.CE_vec=CE_vec;
                return;
            end
            
            % CHANGE
            
            % borrower payoff risk next period
             % Quadrature way
            [nodes, weights] = GaussHermite(11);
            sig2_om_next=sig2_om_next(:);
            rho_om_next=rho_om*ones(size(sig2_om_next));
            % cap on aggregate indexation
            pgr=p_next/p0;
            pgr(pgr>pgr_max)=pgr_max;
            pgr(pgr<pgr_min)=pgr_min;      
            if indexTail == 1
                pgr = pgr + (1-pgr_max);
            end            
            zeta_next=(1 - iota_p + iota_p * pgr) ./ mu_G;  % note: dividing by p0 here
            
            AB_next=nextst(:,2);
            MB_next=nextst(:,3);
            KREO_next = nextst(:, 4);
            KB_next=KBbar-nextst(:,4);
            WI_next=nextst(:, 5);
            BG_next=nextst(:,6);
            QK_next = (1-nuK - (1-ZR_next).*phiK.*lamLTVplus_next).*p_next.*KB_next;
            QA_next = ( (1-tau)*AB_next + delta*(1-ZR_next).*(exp_OmA_B_next.*AB_next) );
            QM_next = ( (1-delta + delta*ZR_next).*MB_next + delta*(1-ZR_next).*(exp_OmM_B_next.*MB_next) );
            

            % Effect on liquidity default rate
            % LIQUIDITY
            QK_liq_next = (1 - nuK) .* liq_thresh .* p_next .* KB_next;
            
            [ZN_next, ZK_next, ZM_next, ZA_next, ...
                Del_NK_next, Del_MK_next, Del_AK_next, Del_KK_next, ...
                Del_NM_next, Del_MM_next, Del_AM_next, Del_KM_next, ...
                Del_NA_next, Del_MA_next, Del_AA_next, Del_KA_next, ...
                ~, ~] ...
                =SAMIntermediaryModel.Mpayoff_logn_quad_liq(QK_next, QA_next, QM_next, QK_liq_next, ...
                KB_next, MB_next, AB_next, p_next, sig2_om_next, rho_om_next, iota_om, zetaombar, ...
                indexAB, indexMB, indexTail,theta, etaB,...
                nodes, weights, false);
           
            % government debt and spending
            bailout_next=F_eps_plus_next - (1-F_eps_next).*(WI_next - zeta_next.*delta.*(1-ZR_next).*(ZA_next.*qA_next.*AB_next+ZM_next.*qM_next.*MB_next));
            TB_next = shareB .* (tauL.*(BG_next + bailout_next));
            
            rent_next=xi_next.*CB_next./((1-xi_next).*KBbar);
            
            Mstar_next=(CB_next - ( (1-tau).*shareB.*Y_next - ZA_next.*(1-tau).*AB_next - ZM_next.*MB_next.*((1-delta) + delta.*ZR_next) ...
                - p_next.*SREO.*KREO_next - rent_next.*KREO_next - nuK.*p_next.*ZK_next.*KB_next - TB_next) )./(ZR_next.*ZN_next);
            
            Kstar_next=(SREO.*KREO_next + ZR_next.*ZK_next.*KB_next)./(ZR_next.*ZN_next);

            
            % New version (DG)
            Del_BC_K_next = internalize_liq * (ZR_next .* (Del_NK_next .* Mstar_next - delta .* Del_MK_next .* MB_next) ...
                - (1.0 - delta) .* Del_MK_next .* MB_next - (1.0 - tau) .* Del_AK_next .* AB_next ...
                - p_next .* (ZR_next .* Del_NK_next .* Kstar_next + (nuK - ZR_next - etaB) .* Del_KK_next .* KB_next));

            Del_BC_M_next = internalize_liq * (ZR_next .* (Del_NM_next .* Mstar_next - delta .* Del_MM_next .* MB_next) ...
                - (1.0 - delta) .* Del_MM_next .* MB_next - (1.0 - tau) .* Del_AM_next .* AB_next ...
                - p_next .* (ZR_next .* Del_NM_next .* Kstar_next + (nuK - ZR_next - etaB) .* Del_KM_next .* KB_next));
            
            Del_BC_A_next = internalize_liq * (ZR_next .* (Del_NA_next .* Mstar_next - delta .* Del_MA_next .* MB_next) ...
                - (1.0 - delta) .* Del_MA_next .* MB_next - (1.0 - tau) .* Del_AA_next .* AB_next ...
                - p_next .* (ZR_next .* Del_NA_next .* Kstar_next + (nuK - ZR_next - etaB) .* Del_KA_next .* KB_next));

            % END LIQUIDITY

            % borrower FOCs
            SDFB=SDF_mat{1};
            indfacAB = (indexAB*zeta_next + 1-indexAB);
            indfacMB = (indexMB*zeta_next + 1-indexMB);
            exp_OmM_B= prnext*( SDFB.*indfacMB/pi.*ZM_next.*(Del_BC_M_next + 1-delta+delta*ZR_next+delta*(1-ZR_next).*exp_OmM_B_next));
            exp_OmA_B= prnext*( SDFB.*indfacAB/pi.*ZA_next.*(Del_BC_A_next + 1-tau+delta*(1-ZR_next).*exp_OmA_B_next));
            exp_OmK_B= prnext*( SDFB.*(rent_next + Del_BC_K_next + p_next.* (etaB.*(1.0-ZK_next) + ZK_next.*(ZR_next-nuK)) + (1-ZR_next).*ZK_next.*exp_OmK_B_next));                        
                         
            % intermediary sector FOC
            SDFI=SDF_mat{2};
            X_next = (1-ZK_next).*KB_next.*(pREO_next-nuREO.*p_next)./MB_next;
            Mpayoff_next=X_next + ZM_next.*(1-delta+delta*ZR_next+delta*(1-ZR_next).*qM_next);
            exp_OmM_I= prnext*( SDFI.*F_eps_next.*indfacMB/pi.*Mpayoff_next);
            Apayoff_next=ZA_next.*(1+delta*(1-ZR_next).*qA_next);
            exp_OmA_I= prnext*( SDFI.*F_eps_next.*indfacAB/pi.*Apayoff_next);
            exp_SDFI=prnext*(SDFI.*F_eps_next)/pi;
			exp_CscrI_next=prnext*(SDFI.*(F_eps_next.*CscrI_next-F_eps_minus_next));
			% no default on REO houses
            exp_pREO_next=prnext*( SDFI.*(rent_next + (SREO-nuREO).*p_next + (1-SREO)*pREO_next));
                                      
            % saver FOC
            SDFS=SDF_mat{3};
            exp_SDFS=prnext*SDFS/pi;
            AB_g = AB_next./indfacAB;
            MB_g = MB_next./indfacMB;
            exp_OmMA_S = prnext*(SDFS.*(indfacMB/pi.*Mpayoff_next + ...
                                         AB_g./MB_g.*indfacAB/pi.*Apayoff_next ));
            
           
            expStruct = struct('exp_OmM_B',exp_OmM_B,...
                   'exp_OmA_B',exp_OmA_B,...
                   'exp_OmK_B',exp_OmK_B,...
                   'exp_OmM_I',exp_OmM_I,...
                   'exp_OmA_I',exp_OmA_I,...
                   'exp_OmMA_S',exp_OmMA_S,...
                   'exp_SDFI',exp_SDFI,...
                   'exp_SDFS',exp_SDFS,...
                   'exp_pREO_next',exp_pREO_next,...
                   'exp_CscrI_next',exp_CscrI_next,...
                   'CE_vec',CE_vec,...
                   'p_next',p_next,...
                   'Apayoff_next',Apayoff_next,...
                   'Mpayoff_next',Mpayoff_next);
               
            Conditional = struct('SDFB',SDFB,...
                   'SDFI',SDFI,...
                   'SDFS',SDFS,...
                   'Mpayoff_next',Mpayoff_next,...
                   'Apayoff_next',Apayoff_next,...
                   'Del_BC_K_next',Del_BC_K_next,...
                   'Del_BC_M_next',Del_BC_M_next);          
               
            expStruct.Conditional = Conditional;  
        end           
        
        
        function [errmat,solmat,condmat,Wshtrans,SDFmat]=calcEEError(obj,pointmat)
            % function to compute Euler equation error at points in state
            % space given by pointmat
            nst=size(obj.Exogenv.pts_all,1);
            
            errmat=zeros(size(pointmat,1),obj.Pfct.Nof);
            solmat=zeros(size(errmat));
            condmat=zeros(size(pointmat,1),obj.NCOND);
            SDFmat=zeros(size(pointmat,1),2*nst);
            Wshtrans=zeros(size(pointmat,1),5*nst);
            
            evaluatePol = @(point)obj.evaluatePol(point);
            calcStateTransition = @(point,soltmp)obj.calcStateTransition(point,soltmp,0);
            calcEquations = @(exst,nextst,soltmp,outstr)obj.calcEquations(exst,nextst,soltmp,outstr,2);    
            % Should be parfor. Use for when debugging only
             parfor i=1:size(errmat,1)
%            for i=1:size(errmat,1)
                point=pointmat(i,:);
                soltmp=evaluatePol(point);
                % transition
                [nextst,outstr]=calcStateTransition(point,soltmp);
                % equations
                [fx,~,V]=calcEquations(point(1),nextst,soltmp,outstr);                                
                qA=exp(soltmp(1));
                qM=exp(soltmp(2));
                qf=exp(soltmp(3));
                p=exp(soltmp(4));
                pREO=exp(soltmp(5));
                Mstar=outstr.addvars.Mstar;
                BSnext=outstr.addvars.BSnext;
                MB_g=outstr.addvars.MB_g;
                AB_g=outstr.addvars.AB_g;                
                normvec=[1,p,1,qA,qM,pREO,qf,qf,p,qM+AB_g/MB_g*qA,qA*AB_g+qM*MB_g,MB_g,BSnext,BSnext];
                condvars=V{1};
                ABtrans=V{2}.AB;                
                MBtrans=V{2}.MB;
                KREOtrans=V{2}.KREO;
                WItrans=V{2}.WI;
                BGtrans=V{2}.BG;
                SDFI=V{3}.SDFI;
                SDFB=V{3}.SDFB;
                errmat(i,:)=fx'./normvec;
                solmat(i,:)=soltmp';
                condmat(i,:)=model.DSGEModel.structToVec(condvars)';
                Wshtrans(i,:) = [ABtrans', MBtrans',KREOtrans',WItrans',BGtrans'];
                SDFmat(i,:) = [SDFI',SDFB'];
            end
            
        end
        
        % simulate model; overwrite method from superclass 
         function [simseries,varnames,errmat,Wshtrans,SDFmat]=simulate(obj,NT,NTini,inistvec,simerror,shmat)
            if length(inistvec)~=obj.Vfct.SSGrid.Ndim
                error('inistvec must be vector of length SSGrid.Ndim');
            end
            
            NTtot=NT+NTini;
            simseries=zeros(NTtot,1+obj.NSTEX+obj.NSTEN+obj.NSOL+obj.NV+obj.NADD);
            
            % if shock matrix wasn't passed in, create it
            preset_path = false;
            if nargin<6
                rng(10,'twister');
                shmat=lhsdesign(NTtot,1,'criterion','correlation'); % realizations of shocks for Markov state
            else 
                preset_path=isinteger(shmat);
            end        
            
            point=inistvec;
                       
            pointmat=zeros(NTtot,length(point));
            
            for t=1:NTtot
               pointmat(t,:)=point; 
               exst=point(1);
               
                % next period's exog. state
               if preset_path
                    exnext=shmat(t);
               else
                    transprob=cumsum(obj.Exogenv.mtrans(exst,:));
                    exnext=find(transprob-shmat(t)>0,1,'first');
               end
               
               % transition to next period
               solvec=obj.evaluatePol(point);
               valvec=obj.evaluateVal(point)';               
               [nextst,outstr]=obj.calcStateTransition(point,solvec,exnext);
%               WIlbvio=nextst(:,5)<obj.Vfct.SSGrid.StateBounds(1,5);
%               nextst(WIlbvio,5)=obj.Vfct.SSGrid.StateBounds(1,5);
%               WIubvio=nextst(:,5)>obj.Vfct.SSGrid.StateBounds(2,5);
%               nextst(WIubvio,5)=obj.Vfct.SSGrid.StateBounds(2,5);
               
               addvec=model.DSGEModel.structToVec(outstr.addvars)';
               % write different categories of variables in one row
               simnext=[point(1),outstr.exstvec',point(2:end),solvec',valvec,addvec];
               if length(simnext)~=size(simseries,2)
                    disp('problem');
               end
               simseries(t,:)=simnext;
               point=nextst;
            end     
            
            simseries=simseries(NTini+1:end,:);
            varnames=[{'exst'}, obj.Ex_names, obj.En_names, obj.Sol_names, obj.V_names, obj.Add_names];
            
            errmat=[];
            Wshtrans=[];%zeros(size(obj.Exogenv.mtrans(simseries(:,1),:)));
            SDFmat=[];%zeros(size(obj.Exogenv.mtrans(simseries(:,1),:)));
            if simerror
                [errmat,~,condmat,Wshtrans,SDFmat]=obj.calcEEError(pointmat);
                errmat=errmat(NTini+1:end,:);
                condmat=condmat(NTini+1:end,:);
                Wshtrans=Wshtrans(NTini+1:end,:);
                SDFmat=SDFmat(NTini+1:end,:);
                simseries=[simseries,condmat];
                varnames=[varnames,obj.Cond_names];
            else
                simseries=[simseries,zeros(NT,obj.NCOND)];                
                varnames=[varnames,strcat(obj.Cond_names,'_nan')];
            end
         end
        
        function [mobj,failedPoints,dist,distT]=polIter(mobj,MAXIT,revisitFailed,printmode,damp,tol_avg)
            gridSt=mobj.Vfct.SSGrid.Pointmat; % use points from BaseGrid here
            NPT=mobj.Vfct.SSGrid.Npt;
            NDIM=mobj.Vfct.SSGrid.Ndim;
            exnpt=size(mobj.Exogenv.pts_perm,1);
            
            % initialize
            resmat=mobj.evaluatePol(gridSt)';
            resmat_prev=resmat;
            
            % split up matrix of points for better output
            gr_points=cell(exnpt,1);
            gr_index=cell(exnpt,2);
            for i=1:exnpt
                grinlog=(gridSt(:,1)==i);
                grind=find(grinlog);
                gr_points{i}=gridSt(grinlog,:);
                gr_index{i,1}=grinlog;
                gr_index{i,2}=grind;
            end
            
            % value function
            VF=mobj.evaluateVal(gridSt)';
            VFnext=zeros(size(VF));
            TF=mobj.evaluateTrans(gridSt)';
            TFnext=zeros(size(VF,1),mobj.Tfct.Nof);
                        
            % control flags
            iter=0;
                        
            disp(' ');
            disp('Starting main loop ...');
            disp(' ');
            while 1
                % counter
                iter=iter+1;
                
                % ===========================================
                % loop over state space
                % ===========================================
                
                % matrix for failed points
                failedPoints=[];
                % transitions
                transmat=mobj.evaluateTrans(gridSt)';
                
                failedPoints_trans_T=zeros(0,size(transmat,2));
                failedPoints_trans_I=[];
                failedPoints_trans_V=zeros(0,size(VF,2));
                
                % a rectangular grid to speed up VF interpolations solutions           
                resmat_startiter = resmat; % for dampening
                % outer loop: all exogenous states
                for ei=1:exnpt
                    tmp_grid=gr_points{ei};
                    tmp_indlog=gr_index{ei,1};
                    tmp_index=gr_index{ei,2};
                    tmp_resmat=resmat(tmp_indlog,:);
                    tmp_resmat_prev=resmat_prev(tmp_indlog,:);
                    
                    tmp_transmat=transmat(tmp_indlog,:);
                    tmp_NPT=size(tmp_index,1); % nb SS pts for exog state ei
                    
                    % Evaluate value functions at transition points
                    transpts=reshape([repmat(1:exnpt,tmp_NPT,1),tmp_transmat],tmp_NPT*exnpt,NDIM);
                    Vtrans = mobj.evaluateVal(transpts)';
                    
                    % index matching for transitions
                    refidx=kron((1:tmp_NPT)',ones(1,exnpt));
                    refidx=refidx(:);                    
                    
                    
                    disp(['State ',num2str(ei)]);

                    [tmp_resmat_new,tmp_VF,tmp_TF,tmp_failed]=mobj.solvePointList(tmp_grid,tmp_resmat,tmp_transmat,...
                        refidx,Vtrans,tmp_resmat_prev,printmode,[]);
                    failedPoints=[failedPoints; tmp_index(tmp_failed)];
                    if revisitFailed
                        failedPoints_trans_T=[failedPoints_trans_T; tmp_transmat(tmp_failed,:)];
                        refidx_failed = ismember(refidx,find(tmp_failed));
                        
                        failedPoints_trans_I=[failedPoints_trans_I; tmp_index(refidx(refidx_failed))];
                        failedPoints_trans_V=[failedPoints_trans_V; Vtrans(refidx_failed,:)];
                    end
                                       
                    resmat_prev(tmp_indlog,:)=tmp_resmat;
                    resmat(tmp_indlog,:)=tmp_resmat_new;
                    VFnext(tmp_indlog,:)=tmp_VF;
                    TFnext(tmp_indlog,:)=tmp_TF;
                end                          
                
                if (revisitFailed && ~isempty(failedPoints))
                    disp( '~~~~~~~~~~~~~~~~~~~');
                    disp(['Revisiting failed points: ',num2str(length(failedPoints)),' add. points ...']);
                    % try to solve at failed points
                    [new_resmat,new_VF,new_TF]=mobj.solvePointListFailed(gridSt,failedPoints,resmat,...
                        failedPoints_trans_T,failedPoints_trans_I,failedPoints_trans_V,...
                        1,printmode,[]);
                    resmat(failedPoints,:)=new_resmat;
                    VFnext(failedPoints,:)=new_VF;
                    TFnext(failedPoints,:)=new_TF;
				end     
                
                
                % approximate functions for next iteration
                mobj=mobj.updateVfct((1-damp)*VFnext+damp*VF);
                mobj=mobj.updateTfct((1-damp)*TFnext+damp*TF);
                
                % convergence criterion (based on points in BaseGrid)
                val_range=4:6;
                VF_val = VF(:,val_range);
                VFnext_val=VFnext(:,val_range);
                [dist,wh]=max(abs(VF_val(:)-VFnext_val(:)));
                [mean_dist,col]=max(abs(mean(VF_val-VFnext_val)));
                [distT,whT]=max(abs(TF(:)-TFnext(:)));                
                [wh_1,wh_2]=ind2sub(size(VFnext_val),wh);
                [whT_1,whT_2]=ind2sub(size(TFnext),whT);
                disp(['-- Iteration: ',num2str(iter),', max distance: ',num2str(dist),' in ',char(mobj.V_names(val_range(1)-1+wh_2)), ...
                    ' at point ',num2str(wh_1),': ',num2str(mobj.Vfct.SSGrid.Pointmat(wh_1,:))]);
                disp(['-- Iteration: ',num2str(iter),', mean distance: ',num2str(mean_dist),' in ',char(mobj.V_names(val_range(1)-1+col))]);
                %disp(num2str(mobj.Vfct.SSGrid.Pointmat(wh_1,:)));
                disp(['-- Iteration: ',num2str(iter),', max T distance: ',num2str(distT),' in col ',num2str(whT_2), ...
                    ' at point ',num2str(whT_1),': ',num2str(mobj.Vfct.SSGrid.Pointmat(whT_1,:))]);
                disp(' ');
%                if dist<0.001 && mean_dist<tol_avg
                if mean_dist<tol_avg
                    disp('Converged.');
                    break;
                elseif iter>=MAXIT
                    disp('Max.iter. exceeded.');
                    break;
                end
                
                % update guess (based on points in BaseGrid)
                VF=VFnext;
                TF=TFnext; 
            end
            
            % resulting policy functions
            mobj=mobj.updatePfct(resmat);       
        end
        
        function [mobj,dist,mean_dist]=convergeVF(mobj,MAXIT,damp,tol_avg)
            
            stateSpace=mobj.Pfct.SSGrid.Pointmat;
            npt=mobj.Pfct.SSGrid.Npt;
            NDIM=mobj.Vfct.SSGrid.Ndim;
            exnpt=size(mobj.Exogenv.pts_perm,1);
                        
            vfindex=4:6;
            tolvf=tol_avg;
            
            iter=0;
            maxiter=MAXIT;
            
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
                newV = oldV;
                nextst=mobj.calcStateTransition(stateSpace(1,:),polmat(1,:),0,transmat(1,:)',oldV(1,:));
                p0=exp(polmat(1,4));
                exptempl = mobj.computeExpectations(stateSpace(1,1),nextst, p0, Vtrans(refidx==1,:));
                
                parfor p=1:npt
                    
                    thispt=stateSpace(p,:);
                    thispol=polmat(p,:);
                    
                    thistrans=transmat(p,:)';
                    vtrans=Vtrans(refidx==p,:);
                    thisvals=oldV(p,:);
                    
                    [nextst,outstr]=mobj.calcStateTransition(thispt,thispol,0,thistrans,thisvals);
                    p0=exp(thispol(4));
                    thisexp=exptempl;
                    thisexp = mobj.computeExpectations(thispt(1),nextst, p0, vtrans, thisexp);
                        
                    
                    %calcEquations(obj,exst,nextst,solvec,instr,mode,varargin)
                    [~,~,thisV]=mobj.calcEquations(thispt(1),nextst,thispol,outstr,3,thisexp);
                    
                    newV(p,vfindex)=thisV{1}(vfindex);
                    
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
                
                mobj=mobj.updateVfct((1-damp)*newV+damp*oldV);
                
                
            end                                    
        end
        
        
        
        function [simseries, varnames] = computeSimulationMoments(obj, simseries, varnames)
                       
            NT_sim = size(simseries,1);
            
            % make HashMap with mapping of names to indices
            indexmap=java.util.HashMap;
            for i=1:length(varnames)
                if isempty(indexmap.get(varnames{i}))
                    indexmap.put(varnames{i},i);
                end
            end
            
            %--------------------------------------------------------------------------
            % transformations of raw model output
            %--------------------------------------------------------------------------

            % list of indices
            loglist=model.HelperCollection.makeListFromNames(indexmap,{'qA','qM','qf','p','cB','cI','cS','pREO','ZR'});
            multlist=model.HelperCollection.makeListFromNames(indexmap,{'lamLTV','lamI'});    

            % conversion of log-values
            simseries(:,loglist)=exp(simseries(:,loglist));
            % conversion of multipliers
            simseries(:,multlist)=max(simseries(:,multlist),0).^(1/3);


            % Function to calculate log growth rates

            G=obj.Params.mu_G;
            loggr = @(var) diff(log(var)) + log(G);            
            params = obj.Params;

            % ---------------------------------------------------------------------
            % state vars (also levels)
            % ---------------------------------------------------------------------
            Y = simseries(:,indexmap.get('Y')); %G = simseries(:,indexmap.get('G'));
            AB = simseries(:,indexmap.get('AB'));
            MB = simseries(:,indexmap.get('MB'));            
            KREO= simseries(:,indexmap.get('KREO'));
            % ---------------------------------------------------------------------
            % borrower debt
            % ---------------------------------------------------------------------
            p= simseries(:,indexmap.get('p'));
            KB= simseries(:,indexmap.get('KB'));
            ZN=simseries(:,indexmap.get('ZN'));
            Drate_ZN= 1-ZN;
            ZK= simseries(:,indexmap.get('ZK'));
            Kstar = simseries(:, indexmap.get('Kstar')); % added by DG
            indexMB=params.indexMB;
            indexAB=params.indexAB;
            indexTail=params.indexTail;
            ZM=simseries(:,indexmap.get('ZM'));
            ZA=simseries(:,indexmap.get('ZA'));
            Drate_ZA=1-ZA;
            Drate_ZM=1-ZM;
            
            % ---------------------------------------------------------------------
            % Maintenance technology
            % ---------------------------------------------------------------------
            
            if params.maint_by_price
                nuK = params.nuK;
                nuREO = params.nuREO;
            else
                nuK = params.nuK_p ./ p;
                nuREO = params.nuREO_p ./ p;
            end
            
            % ---------------------------------------------------------------------
            % wealth distribution
            % ---------------------------------------------------------------------

            WBm= p.*KB - MB;
            pREO= simseries(:,indexmap.get('pREO'));
            WI= simseries(:,indexmap.get('WI'));
            WIm= WI + pREO.*KREO;
            WSm= MB - WIm;

            % ---------------------------------------------------------------------
            % interest rates and returns
            % ---------------------------------------------------------------------

            bind_lamLTV= simseries(:,indexmap.get('lamLTV'))>0;
            bind_lamI= simseries(:,indexmap.get('lamI'))>0;
            delta=params.delta;
            iota_p=params.iota_p;
            indexTail=params.indexTail;
            pgr_max=params.pgr_max;
            pgr_min=params.pgr_min;
            pi=params.pi;
            pgr=p(2:end)./p(1:end-1);
            pgr(pgr>pgr_max)=pgr_max;
            pgr(pgr<pgr_min)=pgr_min;            
            if indexTail == 1
                pgr = pgr + (1-pgr_max);
            end
            zeta = (1 - iota_p + iota_p * pgr) / G;
            rent=simseries(:,indexmap.get('rent'));            
            retp=p(2:end)./p(1:end-1)-1;
            log_p = log(p);
            ZR= simseries(:,indexmap.get('ZR'));
            qM= simseries(:,indexmap.get('qM'));
            qA= simseries(:,indexmap.get('qA'));
            qf= simseries(:,indexmap.get('qf'));
            X= simseries(:,indexmap.get('X'));
            rstar= (1-qM)./qA;
            rbar = AB./MB;
            rD_nom= 1./qf-1;
            rD_real= (1/pi)./qf-1;
            expRP_A= simseries(:,indexmap.get('expRP_A'));
            expRP_M= simseries(:,indexmap.get('expRP_M'));
            expRP_T= simseries(:,indexmap.get('expRP_T'));
            if ~isempty(expRP_A)
                expERP_A=expRP_A - rD_real -1;
                expERP_M=expRP_M - rD_real -1;
                expERP_T=expRP_T - rD_real -1;
            else
                expERP_A=nan(NT_sim,1);
                expERP_M=nan(NT_sim,1);
                expERP_T=nan(NT_sim,1);
            end               
            Lspr_star= rstar - rD_nom;
            Lspr_bar=rbar - rD_nom;
            
            % borrower debt and leverage
            MB_g=simseries(:,indexmap.get('MB_g'));
            AB_g=simseries(:,indexmap.get('AB_g'));            
            MBpayoff_st= (X + ZM.*( 1-delta + delta*ZR + delta*(1-ZR).*qM ));
            MBpayoff=(indexMB*zeta + 1-indexMB)/pi.*MBpayoff_st(2:end);
            ABpayoff_st= ZA.*( 1 + delta*(1-ZR).*qA );
            ABpayoff= (indexAB*zeta + 1-indexAB)/pi.*ABpayoff_st(2:end);
            MBret=MBpayoff./qM(1:end-1) - rD_real(1:end-1) -1;
            ABret=ABpayoff./qA(1:end-1) - rD_real(1:end-1) -1;
            Lret= ( MBpayoff + rstar(1:end-1).*ABpayoff ) - rD_real(1:end-1) -1; % divided by rstar*qA+qM=1
            % total defaulted debt per face value
            deftot=Drate_ZM.* (1- delta + delta*ZR + delta*(1-ZR).*qM) + Drate_ZA .* AB./MB.*(1 + delta*(1-ZR).*qA);
            % total debt value per face value
            debttot=(1- delta + delta*ZR + delta*(1-ZR).*qM) + AB./MB.*(1 + delta*(1-ZR).*qA);
            LGD_ZA = 1 - X./ deftot; 
            % weighted average default rate
            Drate_ZA_avg= deftot./debttot;
            Lrate=Drate_ZA_avg.*LGD_ZA; 
            % without asymmetric indexation
            LGD = 1 - X./ (Drate_ZN.* (1- delta + delta*ZR + delta*(1-ZR).*qM + AB./MB.*(1 + delta*(1-ZR).*qA))); 
            Lrate_ZN=Drate_ZN.*LGD;
            % liquidity defaults
            theta = params.theta;
            Drate_ZA_liq = 1 - simseries(:,indexmap.get('ZA_liq'));
            Drate_ZM_liq = 1 - simseries(:,indexmap.get('ZM_liq'));
            ZK_liq = simseries(:,indexmap.get('ZK_liq'));
            fracliq_ZA = theta * Drate_ZA_liq./Drate_ZA;
            fracliq_ZM = theta * Drate_ZM_liq./Drate_ZM;
            X_liq = X.*(1-ZK_liq)./(1-ZK);
            deftot_liq=Drate_ZM_liq.* (1- delta + delta*ZR + delta*(1-ZR).*qM) + Drate_ZA_liq .* AB./MB.*(1 + delta*(1-ZR).*qA);
            LGD_ZA_liq = 1 - X_liq./ deftot_liq;
            
            % intermediary leverage
            MS_g=simseries(:,indexmap.get('MS_g'));
            fracMS = MS_g./MB_g;
            AS_g = AB_g.*fracMS;
            MI_g = MB_g - MS_g;
            AI_g = AB_g - AS_g;
            AInext = AI_g(1:end-1) .* (indexAB*zeta + 1-indexAB) / pi;
            MInext = MI_g(1:end-1) .* (indexMB*zeta + 1-indexMB) / pi;          

            BInext=simseries(:,indexmap.get('BInext'));
            KREO_g=simseries(:,indexmap.get('KREO_g'));            
            Ifinlvg= BInext./(qM.*MI_g + qA.*AI_g);
            Itotlvg= BInext./(qM.*MI_g + qA.*AI_g + pREO.*KREO_g);           
            
            F_eps=simseries(:,indexmap.get('F_eps'));
            F_eps_minus=simseries(:,indexmap.get('F_eps_minus'));
            JI = qA.*AI_g + qM.*MI_g - qf.*BInext;
            DivI = F_eps.*WI - F_eps_minus - JI;
            Drate_I = 1 - F_eps;            
            DivI_noBR=WI-JI;
            DWL_I = params.zeta*delta*Drate_I(2:end).*(1-ZR(2:end)).*(ZM(2:end).*qM(2:end).*MInext + ZA(2:end).*qA(2:end).*AInext);
            DWL_I = [mean(DWL_I);DWL_I];
            DivREO = (params.SREO-nuREO).*p.*KREO + rent.*KREO - (1-ZK).*pREO.*KB;
            REOinc_c = (params.SREO-nuREO).*p.*KREO + rent.*KREO - nuREO.*(1-ZK).*KB.*p;
            REOexp = (1-ZK).*pREO.*KB;
            maintI = nuK.*params.HI.*p;
            if params.maint_by_price
                retREO = (rent(2:end) + (params.SREO-nuREO).*p(2:end) + (1-params.SREO)*pREO(2:end))./pREO(1:end-1);
            else
                retREO = (rent(2:end) + (params.SREO-nuREO(2:end)).*p(2:end) + (1-params.SREO)*pREO(2:end))./pREO(1:end-1);                
            end
            
            Lstar=simseries(:,indexmap.get('Lstar'));
            WI_noBI = (X(2:end) + ZM(2:end).*(1-delta + delta*ZR(2:end))).*MInext + ZA(2:end).*AInext + delta*(1-ZR(2:end)).*(ZM(2:end).*qM(2:end).*MInext+ZA(2:end).*qA(2:end).*AInext);
            BIbypi = WI_noBI - WI(2:end);
            BIcorrect = BIbypi - BInext(1:end-1)/(G*params.pi);
            netintinc=(X(2:end) +  ZM(2:end).*(1-delta +  delta*ZR(2:end))).*MInext + ZA(2:end).*AInext + qf(2:end).*BInext(2:end) ...
                        - BInext(1:end-1)/(G*params.pi) - Lstar(2:end);
            netintinc_c = netintinc - BIcorrect;        
            
            % ---------------------------------------------------------------------
            % Borrower debt terms
            % ---------------------------------------------------------------------
            
            Mstar=simseries(:,indexmap.get('Mstar'));
            ZN=simseries(:,indexmap.get('ZN'));            
            payB=(1-delta)*ZM.*MB + (1-params.tau)*ZA.*AB;
            issB=ZR.*(ZN.*Mstar - delta*ZM.*MB);
            netisspayB=issB-payB;
            taxB=params.chiB*DWL_I;
            bkLTV= MB./(p.*KB);               
            mkLTV=(qM.*MB+qA.*AB)./(p.*KB);  
            VTI=(p.*KB)./(params.shareB*Y);
            
            % Borr consumption terms (DG)
            net_owned_B = p .* (ZR .* ZN .* Kstar + (nuK - ZR) .* ZK .* KB);
            net_rent_B = rent .* (params.KBbar - KB);
            net_housing_B = net_owned_B + net_rent_B;
            
            % Default terms (DG)
            exp_OmM_B = simseries(:, indexmap.get('exp_OmM_B'));
            exp_OmA_B = simseries(:, indexmap.get('exp_OmA_B'));
            lamLTV = simseries(:, indexmap.get('lamLTV'));
            
            cont_cost_debt = params.delta * (1.0 - ZR) .* (exp_OmM_B .* MB + exp_OmA_B .* AB);
            cont_val_housing = (1.0 - nuK - (1.0 - ZR) .* lamLTV * params.phiK) .* p .* KB;
            default_pay_incentive = (params.delta * ZR + (1.0 - params.delta)) .* MB + (1.0 - params.tau) * AB;
            
            % ---------------------------------------------------------------------
            % Consumption and welfare
            % ---------------------------------------------------------------------
            
            maintKB=nuK.*ZK.*p;
            maintREO=nuREO.*p.*((1-ZK).*KB + KREO);
            maint=maintKB+maintREO;
            maint_byY=maint./Y;
                       
            
            cB=simseries(:,indexmap.get('cB'));
            cI=simseries(:,indexmap.get('cI'));
            cS=simseries(:,indexmap.get('cS'));
            
            C = cB + cI + cS;
            
            VB=simseries(:,indexmap.get('VB'));
            VI=simseries(:,indexmap.get('VIHH'));
            VS=simseries(:,indexmap.get('VS'));
            welfare = VB + VI + VS;

            % Risk sharing
            VCB = VB ./ cB;
            VCI = VI ./ cI;
            VCS = VS ./ cS;
            
            MUB = (1 - params.betaB) * VCB;
            MUI = (1 - params.betaI) * VCI;
            MUS = (1 - params.betaS) * VCS;

            MUBI = log(MUB ./ MUI);
            MUIS = log(MUI ./ MUS);
            MUBS = log(MUB ./ MUS);    
            
            % Compute output ratios
            Y_denom = Y;
            byVar_table = table(KB, KREO, cB, cI, cS, C);
            byVar_varnames=byVar_table.Properties.VariableNames;
            byVar_table = varfun(@(var)var./Y_denom,byVar_table);
            byVar_table.Properties.VariableNames=strcat(byVar_varnames,'_byY');       
                                  
            
            % Growth rates
            grVar_table = table(Y,C,MB,AB,mkLTV,cB,cI,cS,p,issB,payB,WI,DWL_I);
            grVar_varnames=grVar_table.Properties.VariableNames;
            grVar_table = varfun(loggr,grVar_table);
            grVar_table.Properties.VariableNames=strcat(grVar_varnames,'_gr');
            

            
            statevec=simseries(:,1);
            simseries=[simseries(:,2:end),... state vars
                                      WBm, WSm, WIm, fracMS, ... wealth distribution
                                      maintKB, maintREO, maint, maint_byY, maintI, DivREO, REOinc_c, DivI_noBR, REOexp, log_p, ... maintenance
                                      bkLTV, mkLTV, VTI, Drate_ZN, Drate_ZA_avg, LGD, Lrate, LGD_ZA, Lrate_ZN, fracliq_ZA, fracliq_ZM, LGD_ZA_liq, ...
                                      Ifinlvg, Itotlvg, DivI, Drate_I, DWL_I, ...
                                      rstar, rbar, rD_nom, rD_real, bind_lamLTV, bind_lamI, expERP_A, expERP_M, expERP_T,...
                                      Lspr_star, Lspr_bar, payB, issB, netisspayB, taxB,...  debt
                                      welfare, MUBI, MUIS, MUBS, ...
                                      net_owned_B, net_rent_B, net_housing_B, ... %DG
                                      cont_cost_debt, cont_val_housing, default_pay_incentive]; % DG
            simseries = [simseries, byVar_table{:,:}];

            simseries=[simseries(2:end,:),MBret,ABret,Lret,retp,retREO,netintinc,netintinc_c,BIcorrect];
            simseries = [simseries,grVar_table{:,:}];

            varnames_add=[{'WBm','WSm','WIm','fracMS',... wealth distribution
                      'maintKB','maintREO','maint', 'maint_byY', 'maintI', 'DivREO', 'REOinc_c', 'DivI_noBR', 'REOexp', 'logp', ...
                      'bkLTV', 'mkLTV', 'VTI', 'Drate_ZN','Drate_ZA', 'LGD', 'Lrate', 'LGD_ZA', 'Lrate_ZN', 'fracliq_ZA', 'fracliq_ZM', 'LGD_ZA_liq', ...
                      'Ifinlvg', 'Itotlvg', 'DivI', 'Drate_I', 'DWL_I', ...
                      'rstar', 'rbar', 'rD_nom', 'rD_real', 'bind_lamLTV','bind_lamI', 'expERP_A', 'expERP_M', 'expERP_T', ...
                      'Lspr_star', 'Lspr_bar', 'payB', 'issB', 'netisspayB', 'taxB', ... 
                      'welfare', 'MUBI', 'MUIS', 'MUBS', ... consumption and welfare
                      'net_owner_B', 'net_rent_B', 'net_housing_B', ... % DG
                      'cont_cost_debt', 'cont_val_housing', 'default_pay_incentive'}, ...  % DG
                      byVar_table.Properties.VariableNames,...
                      {'MBret','ABret','Lret','retp','retREO','netintinc','netintinc_c','BIcorrect'},...
                      grVar_table.Properties.VariableNames];
            varnames=[varnames(2:end), varnames_add];
            
        end
        
        
    end %of object methods
        
    
    %==============================================================================
    methods (Static)
        % static class-specific methods
        
        function [ZN,ZK,fom,Zbar]=Mpayoff_logn(ombar,mu,sig2,zetaombar)
            lombar=log(ombar);
            sig=sqrt(sig2);
            ZN=1-normcdf((lombar-mu)./sig);
            ZK=normcdf((mu+sig2-lombar)./sig);
            fom=normpdf(lombar,mu,sig);
            Zbar=normcdf((mu+sig2-log(zetaombar))./sig);
        end
  
        
        function XAX = quad_form(X, A)
            XA = mtimesx(X, 't', A); % can use mex-compiled mtimesx here
            XAX = squeeze(mtimesx(XA, X));
        end

        
        % New function including liquidity defaults
        function [ZN, ZK, ZM, ZA, ...
                Del_NK, Del_MK, Del_AK, Del_KK, ...
                Del_NM, Del_MM, Del_AM, Del_KM, ...
                Del_NA, Del_MA, Del_AA, Del_KA, ...
                Z_vals, jacout] ...
                = Mpayoff_logn_quad_liq(QK, QA, QM, QK_liq, ...
                KB, MB, AB, p, sig2, rho, iota_om, zetaombar, ...
                indexAB, indexMB, indexTail,theta, etaB, ...
                nodes, weights, jac)
            
            if (indexTail > 0) || (zetaombar == Inf)
                nnd=length(nodes);
                jacout=[];
                
                sig2_L = iota_om * sig2;
                mu_L = -0.5 * sig2_L;
                sig_L = sqrt(sig2_L);
                
                sig2_U = (1.0 - iota_om) * sig2;
                mu_U = repmat(-0.5 * sig2_U,1,nnd);
                sig_U = repmat(sqrt(sig2_U),1,nnd);
                
                log_om_L = sqrt(2) * sig_L * nodes' + mu_L;
                om_L = exp(log_om_L);
                
                if indexTail == 1
                    om_L_indexed = min(om_L, zetaombar) + (1.0 - zetaombar);
                elseif indexTail == 2
                    lam_L = MB ./ (om_L * p * KB);
                    om_L_indexed = min(1.0, zetaombar ./ lam_L);
                else
                    om_L_indexed = om_L;
                end
                
                om_L_AB = (1.0 - indexAB) + indexAB * om_L_indexed;
                om_L_MB = (1.0 - indexMB) + indexMB * om_L_indexed;
                
                om_U_bar_str = (om_L_MB .* repmat(QM,1,nnd) + om_L_AB .* repmat(QA,1,nnd)) ./ (om_L .* repmat((1.0 + etaB) * QK,1,nnd));
                om_U_bar_str(om_U_bar_str<0)=0;
                
                [ZN_str, ZK_str, ZA_str, ZM_str, jacout_str] ...
                    = SAMIntermediaryModel.computeZ(om_U_bar_str, om_L, om_L_AB, om_L_MB, ...
                        mu_U, sig_U, sig2_U, weights, jac);
                    
                % Liquidity version
                % definition of QK_liq: QK_liq = (1 - nuK) * liq_thresh * p * KB;
                om_U_bar_liq = (om_L_AB .* AB + om_L_MB .* MB) ./ (om_L .* QK_liq);
                om_U_bar_liq(om_U_bar_liq<0)=0;
                
                [ZN_liq, ZK_liq, ZA_liq, ZM_liq, jacout_liq] ...
                    = SAMIntermediaryModel.computeZ(om_U_bar_liq, om_L, om_L_AB, om_L_MB, ...
                        mu_U, sig_U, sig2_U, weights, jac);
                    
                [ZN, ZK, ZA, ZM] = SAMIntermediaryModel.mixZ(theta, ...
                    ZN_liq, ZK_liq, ZA_liq, ZM_liq, ...
                    ZN_str, ZK_str, ZA_str, ZM_str);
                
                Z_vals = [ZN_str, ZK_str, ZM_str, ZA_str, ZN_liq, ZK_liq, ZM_liq, ZA_liq];
%                 frac_liq = (1.0 - ZN_liq) * theta / (1.0 - ZN);
                
                if jac
                    jacout = [jacout_str, jacout_liq];
                end
                
                % omega derivatives
                dom_K = -om_U_bar_liq ./ KB;
                dom_M = om_L_MB ./ (om_L .* QK_liq);
                dom_A = om_L_AB ./ (om_L .* QK_liq);
                
                
                % Z derivatives
                f_om_bar = lognpdf(om_U_bar_liq, mu_U, sig_U);
                dZN = -f_om_bar;
                dZK = -om_L .* om_U_bar_liq .* f_om_bar;
                                
                % Compute \mathcal{I} terms
                I_A = indexAB * om_L + (1.0 - indexAB);
                I_M = indexMB * om_L + (1.0 - indexMB);
                
                dZA = dZN .* I_A; 
                dZM = dZN .* I_M;
                
                
                % Delta terms
                Del_NK = theta * ((dZN .* dom_K) * weights);
                Del_MK = theta * ((dZM .* dom_K) * weights);
                Del_AK = theta * ((dZA .* dom_K) * weights);
                Del_KK = theta * ((dZK .* dom_K) * weights);
                
                Del_NM = theta * ((dZN .* dom_M) * weights);
                Del_MM = theta * ((dZM .* dom_M) * weights);
                Del_AM = theta * ((dZA .* dom_M) * weights);
                Del_KM = theta * ((dZK .* dom_M) * weights);

                Del_NA = theta * ((dZN .* dom_A) * weights);
                Del_MA = theta * ((dZM .* dom_A) * weights);
                Del_AA = theta * ((dZA .* dom_A) * weights);
                Del_KA = theta * ((dZK .* dom_A) * weights);

                
            else
                
                np=length(QK);
                nnd=length(nodes);
                
                ZNcdf_str=zeros(nnd,nnd,np);
                ZKcdf_str=zeros(nnd,nnd,np);
                ZAcdf_str=zeros(nnd,nnd,np);
                ZMcdf_str=zeros(nnd,nnd,np);
                ZNcdf_liq=zeros(nnd,nnd,np);
                ZKcdf_liq=zeros(nnd,nnd,np);
                ZAcdf_liq=zeros(nnd,nnd,np);
                ZMcdf_liq=zeros(nnd,nnd,np);
                
                if jac
                    ZNpdf_str=zeros(nnd,nnd,np);
                    ZKpdf_str=zeros(nnd,nnd,np);
                    ZApdf_str=zeros(nnd,nnd,np);
                    ZMpdf_str=zeros(nnd,nnd,np);
                    ZNpdf_liq=zeros(nnd,nnd,np);
                    ZKpdf_liq=zeros(nnd,nnd,np);
                    ZApdf_liq=zeros(nnd,nnd,np);
                    ZMpdf_liq=zeros(nnd,nnd,np);
                end
                
                dom_K = zeros(nnd, nnd, np);
                dom_M = zeros(nnd, nnd, np);
                dom_A = zeros(nnd, nnd, np);
                dZN = zeros(nnd, nnd, np);
                dZK = zeros(nnd, nnd, np);
                dZA = zeros(nnd, nnd, np);
                dZM = zeros(nnd, nnd, np);
                
                jacout=[];
                for n=1:np
                    % Get R and I components
                    sig2_L = iota_om * sig2(n);
                    mu_L = -0.5 * sig2_L;
                    sig_L = sqrt(sig2_L);
                    sig_L_e = sqrt((1.0 - rho(n)^2) .* sig2_L);
                    
                    sig2_U = (1.0 - iota_om) * sig2(n);
                    mu_U = -0.5 * sig2_U;
                    sig_U = sqrt(sig2_U);
                    
                    % Separate old value and shock
                    log_om_L_lag = sqrt(2) * sig_L * nodes + mu_L;
                    e_L = sqrt(2) * sig_L_e * nodes;
                    
                    % (old x shock)
                    log_om_L_lag_mat = repmat(log_om_L_lag, 1, nnd);
                    e_L_mat = repmat(e_L', nnd, 1);
                    log_om_L = (1.0 - rho(n)) * mu_L + rho(n) * log_om_L_lag_mat + e_L_mat;
                    
                    om_L = exp(log_om_L);
                    om_L_lag = exp(log_om_L_lag_mat);
                    om_L_hat = min(om_L, zetaombar * om_L_lag);
                    
                    om_L_hat_AB = (1.0 - indexAB) + indexAB * om_L_hat;
                    om_L_hat_MB = (1.0 - indexMB) + indexMB * om_L_hat;
                    
                    om_U_bar_str = (om_L_hat_MB * QM(n) + om_L_hat_AB * QA(n)) ./ ((1.0 + etaB) * QK(n) * om_L);
                    om_U_bar_str(om_U_bar_str<0)=0;
                    log_om_U_bar_str = log(om_U_bar_str);
                    
                    % Strategic version
                    [ZNcdf_str(:, :, n), ZKcdf_str(:, :, n), ZAcdf_str(:, :, n), ZMcdf_str(:, :, n)] = ...
                        SAMIntermediaryModel.computeZcdfs(log_om_U_bar_str, om_L_hat_AB, om_L_hat_MB, mu_U, sig2_U, sig_U);
                        
                    if jac
                        [ZNpdf_str(:, :, n), ZKpdf_str(:, :, n), ZApdf_str(:, :, n), ZMpdf_str(:, :, n)] = ...
                            SAMIntermediaryModel.computeZpdfs(log_om_U_bar_str, om_L_hat_AB, om_L_hat_MB, mu_U, sig2_U, sig_U);
                    end
                    
                    % Liquidity version
                    om_U_bar_liq = (om_L_hat_AB .* AB(n) + om_L_hat_MB .* MB(n)) ./ (QK_liq(n) .* om_L);
                    om_U_bar_liq(om_U_bar_liq<0)=0;
                    log_om_U_bar_liq = log(om_U_bar_liq);
                    
                    [ZNcdf_liq(:, :, n), ZKcdf_liq(:, :, n), ZAcdf_liq(:, :, n), ZMcdf_liq(:, :, n)] = ...
                        SAMIntermediaryModel.computeZcdfs(log_om_U_bar_liq, om_L_hat_AB, om_L_hat_MB, mu_U, sig2_U, sig_U);
                        
                    if jac
                        [ZNpdf_liq(:, :, n), ZKpdf_liq(:, :, n), ZApdf_liq(:, :, n), ZMpdf_liq(:, :, n)] = ...
                            SAMIntermediaryModel.computeZpdfs(log_om_U_bar_liq, om_L_hat_AB, om_L_hat_MB, mu_U, sig2_U, sig_U);
                    end
                    

                    dom_K(:, :, n) = -om_U_bar_liq ./ KB(n);
                    dom_M(:, :, n) = om_L_hat_MB ./ (QK_liq(n) .* om_L);
                    dom_A(:, :, n) = om_L_hat_AB ./ (QK_liq(n) .* om_L);
                    
                    
                    f_om_bar = lognpdf(om_U_bar_liq, mu_U, sig_U);
                    dZN(:, :, n) = -f_om_bar;
                    dZK(:, :, n) = -om_L .* om_U_bar_liq .* f_om_bar;
                    
                    % Compute \mathcal{I} terms
                    I_A = indexAB * om_L + (1.0 - indexAB);
                    I_M = indexMB * om_L + (1.0 - indexMB);

                    dZA(:, :, n) = dZN(:, :, n) .* I_A;  
                    dZM(:, :, n) = dZN(:, :, n) .* I_M;

                    
                end
                % Compute fractions of non-defaulters and retained shares
                ZN_str = SAMIntermediaryModel.quad_form(weights, ZNcdf_str);
                ZK_str = SAMIntermediaryModel.quad_form(weights, ZKcdf_str);
                ZA_str = SAMIntermediaryModel.quad_form(weights, ZAcdf_str);
                ZM_str = SAMIntermediaryModel.quad_form(weights, ZMcdf_str);
                
                ZN_liq = SAMIntermediaryModel.quad_form(weights, ZNcdf_liq);
                ZK_liq = SAMIntermediaryModel.quad_form(weights, ZKcdf_liq);
                ZA_liq = SAMIntermediaryModel.quad_form(weights, ZAcdf_liq);
                ZM_liq = SAMIntermediaryModel.quad_form(weights, ZMcdf_liq);
                
                [ZN, ZK, ZA, ZM] = SAMIntermediaryModel.mixZ(theta, ...
                    ZN_liq, ZK_liq, ZA_liq, ZM_liq, ...
                    ZN_str, ZK_str, ZA_str, ZM_str);
                
                Z_vals = [ZN_str, ZK_str, ZM_str, ZA_str, ZN_liq, ZK_liq, ZM_liq, ZA_liq];
                
                if jac
                    ZNp_str = SAMIntermediaryModel.quad_form(weights, ZNpdf_str);
                    ZKp_str = SAMIntermediaryModel.quad_form(weights, ZKpdf_str);
                    ZAp_str = SAMIntermediaryModel.quad_form(weights, ZApdf_str);
                    ZMp_str = SAMIntermediaryModel.quad_form(weights, ZMpdf_str);
                    jacout_str=[ZNp_str, ZKp_str, ZAp_str, ZMp_str];
                    
                    ZNp_liq = SAMIntermediaryModel.quad_form(weights, ZNpdf_liq);
                    ZKp_liq = SAMIntermediaryModel.quad_form(weights, ZKpdf_liq);
                    ZAp_liq = SAMIntermediaryModel.quad_form(weights, ZApdf_liq);
                    ZMp_liq = SAMIntermediaryModel.quad_form(weights, ZMpdf_liq);
                    jacout_liq=[ZNp_liq, ZKp_liq, ZAp_liq, ZMp_liq];
                    
                    jacout = [jacout_str, jacout_liq];
                end
                
                Del_NK = theta * SAMIntermediaryModel.quad_form(weights,  dZN .* dom_K);
                Del_KK = theta * SAMIntermediaryModel.quad_form(weights,  dZK .* dom_K);
                Del_AK = theta * SAMIntermediaryModel.quad_form(weights,  dZA .* dom_K);
                Del_MK = theta * SAMIntermediaryModel.quad_form(weights,  dZM .* dom_K);
                
                Del_NM = theta * SAMIntermediaryModel.quad_form(weights,  dZN .* dom_M);
                Del_KM = theta * SAMIntermediaryModel.quad_form(weights,  dZK .* dom_M);
                Del_AM = theta * SAMIntermediaryModel.quad_form(weights,  dZA .* dom_M);
                Del_MM = theta * SAMIntermediaryModel.quad_form(weights,  dZM .* dom_M);

                Del_NA= theta * SAMIntermediaryModel.quad_form(weights,  dZN .* dom_A);
                Del_KA = theta * SAMIntermediaryModel.quad_form(weights,  dZK .* dom_A);
                Del_AA = theta * SAMIntermediaryModel.quad_form(weights,  dZA .* dom_A);
                Del_MA = theta * SAMIntermediaryModel.quad_form(weights,  dZM .* dom_A);
                
            end
            
        end
        
        
        
        function [ZN, ZK, ZA, ZM, jacout] ...
                = computeZ(om_U_bar, om_L, om_L_AB, om_L_MB, mu_U, sig_U, sig2_U, weights, jac)
            
            % Compute normal CDF terms
            ZNcdf =normcdf( (log(om_U_bar) - mu_U) ./ sig_U );
            ZKcdf = normcdf((mu_U + sig2_U - log(om_U_bar)) ./ sig_U );
            
            % Compute fractions of non-defaulters and retained shares
            ZN_conditional = (1.0 - ZNcdf);
            ZN = (1.0 - ZNcdf) * weights;
            ZK = (ZKcdf .* om_L) * weights;
            ZA = (ZN_conditional .* om_L_AB) * weights;
            ZM = (ZN_conditional .* om_L_MB) * weights;
            
            if jac
                ZNpdf = normpdf( (log(om_U_bar) - mu_U) ./ sig_U ) ./ sig_U;
                ZKpdf = -normpdf((mu_U + sig2_U - log(om_U_bar)) ./ sig_U ) ./ sig_U;
                ZNp_conditional = - ZNpdf;
                ZNp = - ZNpdf * weights;
                ZKp = (ZKpdf .* om_L) * weights;
                ZAp = (ZNp_conditional .* om_L_AB) * weights;
                ZMp = (ZNp_conditional .* om_L_MB) * weights;
                jacout = [ZNp, ZKp, ZAp, ZMp];
            else
                jacout = [];
            end
            
        end
        
        function [ZNcdf, ZKcdf, ZAcdf, ZMcdf] ...
                = computeZcdfs(log_om_U_bar, om_L_hat_AB, om_L_hat_MB, mu_U, sig2_U, sig_U)   
            
            ZNcdf=1 - normcdf((log_om_U_bar - mu_U) ./ sig_U);
            ZKcdf=normcdf((mu_U + sig2_U - log_om_U_bar) ./ sig_U);
            ZAcdf=ZNcdf.* om_L_hat_AB;
            ZMcdf=ZNcdf.* om_L_hat_MB;
            
        end
            
        function [ZNpdf, ZKpdf, ZApdf, ZMpdf] ...
                = computeZpdfs(log_om_U_bar, om_L_hat_AB, om_L_hat_MB, mu_U, sig2_U, sig_U)   
                
            ZNpdf=-normpdf((log_om_U_bar - mu_U) ./ sig_U) ./ sig_U;
            ZKpdf=-normpdf((mu_U + sig2_U - log_om_U_bar) ./ sig_U) ./ sig_U;
            ZApdf=ZNpdf.* om_L_hat_AB;
            ZMpdf=ZNpdf.* om_L_hat_MB;
            
        end
        
        function [ZN, ZK, ZA, ZM] = mixZ(theta, ...
                ZN_liq, ZK_liq, ZA_liq, ZM_liq, ...
                ZN_str, ZK_str, ZA_str, ZM_str)
            
            ZN = theta * ZN_liq + (1.0 - theta) * ZN_str;
            ZK = theta * ZK_liq + (1.0 - theta) * ZK_str;
            ZA = theta * ZA_liq + (1.0 - theta) * ZA_str;
            ZM = theta * ZM_liq + (1.0 - theta) * ZM_str;
                    
        end
        
        function [ZR,fkappa]=Refi(kappabar,mu,sig)
            ZR=1./(4*(1+exp(-(kappabar-mu)./sig)));
            fkappa=exp(-(kappabar-mu)./sig)./(sig.*(1+exp(-(kappabar-mu)./sig)).^2);
        end             
        
        function [J,fxout] = constructJacobian(params)
            syms Y sig2_om xi
            syms AB MB KREO WI BG
            
            syms qA qM qf p pREO cB cI cS ZR MS_g
            syms lamLTVplus lamLTVminus lamIplus lamIminus muSplus muSminus lamSplus lamSminus
            
            syms exp_OmA_B exp_OmM_B exp_OmK_B exp_OmA_I exp_OmM_I exp_OmMA_S exp_SDFI exp_SDFS exp_pREO_next 
            syms OmA OmM CscrI
            syms nodes
            
            % unpack relevant parameters
            shareB=params.shareB;
            shareI=params.shareI;
            shareS=params.shareS;
            nuK=params.nuK;
            nuREO=params.nuREO;
            KBbar=params.KBbar;
            SREO=params.SREO;
            iota_om=params.iota_om;
            psi=params.psi;
            xiL=params.xiL;
            phiK=params.phiK;
            phiI=params.phiI;
            tau=params.tau;
            delta=params.delta;
            HI=params.HI;           
            HS=params.HS;           
            mu_kappa=params.mu_kappa;
            sig_kappa=params.sig_kappa;
            sigma_eps=params.sigma_eps;
            zetaDWL=params.zeta;
            indexMB=params.indexMB;
            indexAB=params.indexAB;
            indexTail=params.indexTail;
            tauL=params.tauL;
            theta=params.theta;
            etaB=params.etaB;
            liq_thresh=params.liq_thresh;
            psiS0=params.psiS0;
            psiS1=params.psiS1;
            
            if params.maint_by_price
                nuK = params.nuK;
                nuREO = params.nuREO;
            else
                nuK = params.nuK_p / p;
                nuREO = params.nuREO_p / p;
            end
            if params.rebate_maint
                nuK_budget = 0;
            else
                nuK_budget = nuK;
            end                
            
            
            syms fZZN(log_omega) fZZK(log_omega) fZZA(log_omega) fZZM(log_omega)
            
            % Get R and I components
            sig2_L = iota_om * sig2_om;
            mu_L = -0.5 * sig2_L;
            sig_L = sqrt(sig2_L);
            
            % Compute omega_U threshold for given omega_L
            log_om_L = sqrt(2) * sig_L * nodes + mu_L;
            om_L = exp(log_om_L);
            om_L_AB = (1.0 - indexAB) + indexAB * om_L;
            om_L_MB = (1.0 - indexMB) + indexMB * om_L;

            KB=KBbar-KREO;
            QK = (1-nuK - (1-ZR)*phiK*lamLTVplus)*p*KB;
            QA = ( (1-tau)*AB + delta*(1-ZR)*(OmA*AB) );
            QM = ( (1-delta + delta*ZR)*MB + delta*(1-ZR)*(OmM*MB) );            
            
            % LIQUIDITY
            QK_liq = (1-nuK) * liq_thresh * p * KB;

            log_om_U_bar_str = log((om_L_MB * QM + om_L_MB * QA) ./ ((1.0 + etaB) * QK * om_L));
            log_om_U_bar_liq = log((om_L_AB * AB + om_L_MB * MB) ./ (QK_liq * om_L));

            ZN_str = fZZN(log_om_U_bar_str);
            ZK_str = fZZK(log_om_U_bar_str);
            ZA_str = fZZA(log_om_U_bar_str);
            ZM_str = fZZM(log_om_U_bar_str);

            ZN_liq = fZZN(log_om_U_bar_liq);
            ZK_liq = fZZK(log_om_U_bar_liq);
            ZA_liq = fZZA(log_om_U_bar_liq);
            ZM_liq = fZZM(log_om_U_bar_liq);

            [ZN, ZK, ZA, ZM] = SAMIntermediaryModel.mixZ(theta, ...
                ZN_liq, ZK_liq, ZA_liq, ZM_liq, ...
                ZN_str, ZK_str, ZA_str, ZM_str);
            % END LIQUIDITY
                        		
            % intermediary default
            syms fF_eps(V) fF_eps_cond(V)			
			VI=WI+CscrI;
            F_eps=fF_eps(VI);
            F_eps_cond=fF_eps_cond(VI);
            F_eps_minus=-sigma_eps*F_eps_cond;
            F_eps_plus=sigma_eps*F_eps_cond;	

            % government debt and spending
			bailout=F_eps_plus - (1-F_eps)*(WI - zetaDWL*delta*(1-ZR)*(ZA*qA*AB+ZM*qM*MB));
            Tvec=[shareB, shareI, shareS]' * (tauL*(BG + bailout));
            BGnext=(1-tauL)*(BG + bailout)/qf; % gov debt EOP

			
            % intraperiod FOC for current interest rate
            rstar=(1-qM)/qA;
            rbar=AB/MB;
            
            % borrower budget            
            rent=xi*cB/((1-xi)*KBbar);
            Mstar=(cB - ( (1-tau)*shareB*Y - ZA*(1-tau)*AB - ZM*(1-delta)*MB - delta*ZR*ZM*MB - p*SREO*KREO - rent*KREO - nuK_budget.*p*ZK*KB - Tvec(1)) )/(ZR*ZN);
            Kstar=(SREO*KREO + ZR*ZK*KB)/(ZR*ZN);
            AB_g=ZR*ZN*rstar*Mstar+delta*(1-ZR)*ZA*AB;
            MB_g=ZR*ZN*Mstar+delta*(1-ZR)*ZM*MB;
                      
            % mortgage market clearing
            AS_g = AB_g*MS_g/MB_g;
            MI_g = MB_g - MS_g;
            AI_g = AB_g - AS_g;
            
            % intermediary and saver budget
            X = (1-ZK)*KB*(pREO-nuREO*p)/MB;
            debtB = (X + ZM*(1-delta + delta*ZR))*MB + ZA*AB + delta*(1-ZR)*(ZM*qM*MB+ZA*qA*AB);
            wealthS = debtB - WI + BG; % total debt held by savers BOP
            BSnext = -(cS - ((1-tau)*shareS*Y + wealthS - qM*MS_g - qA*AS_g - nuK_budget.*p*HS - Tvec(3)))/qf;  % saver debt EOP
            cy = xiL*cS./((1-xi-xiL).*BSnext);
            % impose MC in mortgage market to compute intermediary budget
            KREO_g = (1-SREO)*KREO + (1-ZK)*KB;
            BInext=(cI - ( F_eps*WI - F_eps_minus + (1-tau)*shareI*Y - qM*MI_g - qA*AI_g ...
                + (SREO-nuREO)*p*KREO + rent*KREO - pREO*(1-ZK)*KB - nuK_budget.*p*HI - Tvec(2)))/qf;
            
                        
            % compute (sym) current U's to normalize back SDFs
            U_vecB=cB^(1-xi)*KBbar^xi;
            U_vecI=cI^(1-xi)*HI^xi;
            U_vecS=cS^(1-xi-xiL)* HS^xi * BSnext^xiL;
            U_normB=U_vecB^(1-1/psi)/cB;  
            U_normI=U_vecI^(1-1/psi)/cI;   
            U_normS=U_vecS^(1-1/psi)/cS;
            U_normB_p=U_normB; % if indexation, also divide expectations by p
            U_normI_p=U_normI;

            
            % borrower FOCs
            exp_OmM_B_norm=exp_OmM_B/U_normB_p;
            exp_OmA_B_norm=exp_OmA_B/U_normB_p;
            exp_OmK_B_norm=exp_OmK_B/U_normB;          
            kappabar = (1-exp_OmM_B_norm-rbar*exp_OmA_B_norm)*(1-delta*ZM*MB/(ZN*Mstar)) ...
                             - exp_OmA_B_norm*(rstar-rbar) - p*phiK*lamLTVplus*(ZN*Kstar-ZK*KB)/(ZN*Mstar);
            fx1 = 1-lamLTVplus - exp_OmM_B_norm - rstar*exp_OmA_B_norm;
            fx2 = p*(1-lamLTVplus*phiK) - exp_OmK_B_norm;
            fx3 = ZR - 1/(4*(1+exp(-(kappabar-mu_kappa)./sig_kappa)));

            
            % lender FOCs
            fx4=qA*(1-phiI*lamIplus) - exp_OmA_I/U_normI_p;
            fx5=qM*(1-phiI*lamIplus) - exp_OmM_I/U_normI_p;
            fx6=pREO - exp_pREO_next/U_normI;
            fx7=qf -lamIplus - exp_SDFI/U_normI;                
            
            % saver FOC
            rbar_g=AB_g/MB_g;
            fx8=qf - cy - exp_SDFS/U_normS - lamSplus;
            fx9=qM + rbar_g*qA + psiS0*(MS_g)^(psiS1-1) - exp_OmMA_S/U_normS - muSplus;
                        
            % borrower constraint
            fx10=phiK*p*Kstar - Mstar - lamLTVminus; 
            
            % intermediary constraint
            fx11=phiI*(1-MS_g/MB_g)*(qA*AB_g+qM*MB_g) - BInext - lamIminus;
            
            % saver no-shorting constraint
            fx12=MS_g - muSminus;      
            fx13=BSnext - lamSminus;
            
            % bond market clearing
            fx14=BInext + BGnext - BSnext;
            
            fx=[fx1;fx2;fx3;fx4;fx5;fx6;fx7;fx8;fx9;fx10;fx11;fx12;fx13;fx14];
            
            Jacobian = jacobian(fx,[qA,qM,qf,p,pREO,cB,cI,cS,ZR,MS_g,lamLTVplus,lamLTVminus,lamIplus,lamIminus,muSplus,muSminus,lamSplus,lamSminus]);
            
            syms ZZN_str ZZK_str ZZA_str ZZM_str dZZN_str dZZK_str dZZA_str dZZM_str ...
                ZZN_liq ZZK_liq ZZA_liq ZZM_liq dZZN_liq dZZK_liq dZZA_liq dZZM_liq ...
                F_eps F_eps_cond dF_eps dF_eps_cond
                        
            
            subs_out = {...
                fZZN(log_om_U_bar_str), fZZK(log_om_U_bar_str), fZZA(log_om_U_bar_str), fZZM(log_om_U_bar_str), ...
                fZZN(log_om_U_bar_liq), fZZK(log_om_U_bar_liq), fZZA(log_om_U_bar_liq), fZZM(log_om_U_bar_liq), ...
                subs(diff(fZZN, log_omega), log_omega, log_om_U_bar_str), ...
                subs(diff(fZZK, log_omega), log_omega, log_om_U_bar_str), ...
                subs(diff(fZZM, log_omega), log_omega, log_om_U_bar_str), ...
                subs(diff(fZZA, log_omega), log_omega, log_om_U_bar_str), ...
                subs(diff(fZZN, log_omega), log_omega, log_om_U_bar_liq), ...
                subs(diff(fZZK, log_omega), log_omega, log_om_U_bar_liq), ...
                subs(diff(fZZM, log_omega), log_omega, log_om_U_bar_liq), ...
                subs(diff(fZZA, log_omega), log_omega, log_om_U_bar_liq), ...
                fF_eps(VI), fF_eps_cond(VI), ...
                subs(diff(fF_eps, V), V, VI), ...
                subs(diff(fF_eps_cond, V), V, VI) ...
                };
            
            subs_in = {...
                ZZN_str, ZZK_str, ZZA_str, ZZM_str, ...
                ZZN_liq, ZZK_liq, ZZA_liq, ZZM_liq, ...
                dZZN_str, dZZK_str, dZZA_str, dZZM_str, ...
                dZZN_liq, dZZK_liq, dZZA_liq, dZZM_liq, ...
                F_eps, F_eps_cond, dF_eps, dF_eps_cond ...
                };
            
            fx = subs(fx, subs_out, subs_in);
            Jacobian = subs(Jacobian, subs_out, subs_in);
            
            % Additional work to get rid of D terms
            syms TEMP_str TEMP_liq;
            Jacobian = subs(Jacobian, log_om_U_bar_str, TEMP_str);
            Jacobian = subs(Jacobian, log_om_U_bar_liq, TEMP_liq);
            for ii = 1 : size(Jacobian, 1)
                tic
                for jj = 1 : size(Jacobian, 2)
                    Jac_ij = sprintf('%s', Jacobian(ii, jj));
                    if contains(Jac_ij, 'D(fZZK)') ...
                            || contains(Jac_ij, 'D(fZZN)') ...
                            || contains(Jac_ij, 'D(fZZA)') ...
                            || contains(Jac_ij, 'D(fZZM)')
                        Jac_ij = strrep(Jac_ij, 'D(fZZK)(TEMP_str)', 'dZZK_str');
                        Jac_ij = strrep(Jac_ij, 'D(fZZN)(TEMP_str)', 'dZZN_str');
                        Jac_ij = strrep(Jac_ij, 'D(fZZM)(TEMP_str)', 'dZZM_str');
                        Jac_ij = strrep(Jac_ij, 'D(fZZA)(TEMP_str)', 'dZZA_str');
                        Jac_ij = strrep(Jac_ij, 'D(fZZK)(TEMP_liq)', 'dZZK_liq');
                        Jac_ij = strrep(Jac_ij, 'D(fZZN)(TEMP_liq)', 'dZZN_liq');
                        Jac_ij = strrep(Jac_ij, 'D(fZZM)(TEMP_liq)', 'dZZM_liq');
                        Jac_ij = strrep(Jac_ij, 'D(fZZA)(TEMP_liq)', 'dZZA_liq');
                        Jacobian(ii, jj) = str2sym(Jac_ij);
                    end
                end
                toc
            end
            Jacobian = subs(Jacobian, TEMP_str, log_om_U_bar_str);
            Jacobian = subs(Jacobian, TEMP_liq, log_om_U_bar_liq);
            
            
            order=[Y, sig2_om, xi, AB, MB, KREO, WI, BG,...
                qA, qM, qf, p, pREO, cB, cI, cS, ZR, MS_g,...
                lamLTVplus,lamLTVminus,lamIplus,lamIminus,muSplus,muSminus,lamSplus,lamSminus,...
                ZZN_str,ZZK_str,ZZA_str,ZZM_str,dZZN_str,dZZK_str,dZZA_str,dZZM_str,...
                ZZN_liq,ZZK_liq,ZZA_liq,ZZM_liq,dZZN_liq,dZZK_liq,dZZA_liq,dZZM_liq,...
                F_eps, F_eps_cond, dF_eps, dF_eps_cond,...				
                exp_OmA_B, exp_OmM_B, exp_OmK_B, exp_OmA_I, exp_OmM_I, exp_OmMA_S, exp_SDFI, exp_SDFS, exp_pREO_next,...
                OmA, OmM, nodes];
            J=matlabFunction(Jacobian,'Vars',order); % convert symbolic function to numeric Matlab function (Jacobian)
            fxout=matlabFunction(fx,'Vars',order); % convert symbolic function to numeric Matlab function (residuals)                       
        end        
                
      function [fx,stvals]=compStSt(sol,params,omstate,print)
            % unpack relevant parameters
            betaB=params.betaB;
            betaI=params.betaI;
            betaS=params.betaS;
            shareB=params.shareB;
            shareI=params.shareI;
            shareS=params.shareS;
            nuREO=params.nuREO;
            KBbar=params.KBbar;
            SREO=params.SREO;
            iota_om=params.iota_om;
            psi=params.psi;
            pi=params.pi;
            xivec=[params.xihi,params.xilo];
            xi=xivec(omstate);
            phiK=params.phiK;
            phiI=params.phiI;
            tau=params.tau;
            delta=params.delta;
            xiL=params.xiL;     
            if xiL==0
                xiL=.001; % just to compute steady state
            end
            mu_kappa=params.mu_kappa;
            sig_kappa=params.sig_kappa;
            sig2_om=log(1+[params.sigom_lo,params.sigom_hi].^2);
            sig2_om=sig2_om(omstate);
            g=params.g;
            mu_G=params.mu_G;
            shareHB=params.shareHB;
            shareHI=params.shareHI;
            shareHS=params.shareHS;
            Y=params.mu_Y;
			sigma_eps=params.sigma_eps;
			zetaDWL=params.zeta;
            zetaombar=params.pgr_max;
            indexMB=params.indexMB;
            indexAB=params.indexAB;
            indexTail=params.indexTail;
            tauL=params.tauL;
            theta=params.theta;
            liq_thresh=params.liq_thresh;
            etaB=params.etaB;
            internalize_liq=params.internalize_liq;
            nuK = params.nuK;
            if isfield(params,'rho_om')
                rho_om=params.rho_om;
            else
                rho_om=0.95;
            end
                      
            % solution vars
            p=exp(sol(1));
            cy=exp(sol(2));
            cB=exp(sol(3));
            kappabar=exp(sol(4));
			VI=exp(sol(5));
            QA=exp(sol(6));
            QM=exp(sol(7));
            QK=exp(sol(8));
            KB=exp(sol(9));
            MB=exp(sol(10));
            AB=exp(sol(11));
            
            % Maintenance technology
            if params.rebate_maint
                nuK_budget = 0;
            else
                nuK_budget = nuK;
            end                

            % refi and default
            ZR=SAMIntermediaryModel.Refi(kappabar,mu_kappa,sig_kappa);
            
            [nodes, weights] = GaussHermite(11);
            % LIQUIDITY DEFAULT VERSION
            QK_liq = (1 - nuK) * liq_thresh * p * KB;
            
            [ZN, ZK, ZM, ZA, ...
                Del_NK, Del_MK, Del_AK, Del_KK, ...
                Del_NM, Del_MM, Del_AM, Del_KM, ...
                Del_NA, Del_MA, Del_AA, Del_KA, ...
                Z_vals, ~] ...
                =SAMIntermediaryModel.Mpayoff_logn_quad_liq(QK, QA, QM, QK_liq, ...
                KB, MB, AB, p, sig2_om, rho_om, iota_om, zetaombar, ...
                indexAB, indexMB, indexTail,theta, etaB,...
                nodes, weights, false);
            
            KREO=KB*(1-ZK)/SREO;
            Kstar=(SREO*KREO+ZR*ZK*KB)/(ZR*ZN);


            % short-term bond prices
            betaB_g=betaB*exp(-g/psi);
            betaI_g=betaI*exp(-g/psi);
            betaS_g=betaS*exp(-g/psi);
            zeta=1/pi;            
            qf=cy + betaS_g*zeta;
            rD=1/qf-1;


			% intermediary default
            F_eps=normcdf(VI,0,sigma_eps);
            betaI_g_F=betaI_g*F_eps;
            F_eps_minus=-sigma_eps*normpdf(VI/sigma_eps);
            F_eps_plus=sigma_eps*normpdf(VI/sigma_eps);           
            lamI=qf-betaI_g_F*zeta;

			
            % mortgages
            Mstar=phiK*p*Kstar;
            Lstar=Mstar*ZR*ZN;
            % LIQUIDITY derivate terms
            Del_BC_K = internalize_liq * (ZR * (Del_NK * Mstar - delta * Del_MK * MB) ...
                - (1.0 - delta) * Del_MK * MB - (1.0 - tau) * Del_AK * AB ...
                - p * (ZR * Del_NK * Kstar + (nuK - ZR - etaB) * Del_KK * KB));
            
            Del_BC_M = internalize_liq * (ZR * (Del_NM * Mstar - delta * Del_MM * MB) ...
                - (1.0 - delta) * Del_MM * MB - (1.0 - tau) * Del_AM * AB ...
                - p * (ZR * Del_NM * Kstar + (nuK - ZR - etaB) * Del_KM * KB));
            
            Del_BC_A = internalize_liq * (ZR * (Del_NA * Mstar - delta * Del_MA * MB) ...
                - (1.0 - delta) * Del_MA * MB - (1.0 - tau) * Del_AA * AB ...
                - p * (ZR * Del_NA * Kstar + (nuK - ZR - etaB) * Del_KA * KB));

            % END LIQUIDITY
            
            rent=xi*cB/((1-xi)*KBbar);
            pREO=betaI_g*(rent+(SREO-nuREO)*p)/(1-betaI_g*(1-SREO));
            X=(1-ZK)*(pREO-nuREO*p)*KB/MB;
            qA=betaI_g_F*zeta*ZA/(1-phiI*lamI-betaI_g_F*zeta*delta*ZA*(1-ZR));
            qM=betaI_g_F*zeta*(X +ZM*(1-delta+delta*ZR))/(1-phiI*lamI-betaI_g_F*zeta*delta*ZM*(1-ZR)); % Switching to intermediary version of Del_BC_M and putting inside ZM term
            rstar=(1-qM)/qA;
            
            OmA=betaB_g*zeta*ZA*(1-tau)/(Del_BC_A + 1-betaB_g*zeta*delta*ZA*(1-ZR));
            OmM=betaB_g*zeta*ZM*(Del_BC_M+1-delta+delta*ZR)/(1-betaB_g*zeta*delta*ZM*(1-ZR)); % Adding Del_BC_M (DG)
            OmK=betaB_g*(rent + ZK*p*(ZR-nuK))/(1-betaB_g*ZK*(1-ZR));
            lamLTV=1-OmM-rstar*OmA;
            rbar=AB/MB;
			
            % intermediary wealth, consumption
			MI=MB;
            AI=AB;
            BI=phiI*(qA*AI+qM*MI);
            WI = (X + ZM*(1-delta + delta*ZR))*MI + ZA*AI + delta*(1-ZR)*(ZM*qM*MI+ZA*qA*AI) - BI/pi;            
			CscrI = -betaI_g*F_eps_minus/(1-betaI_g_F);

            
            % gov spending, bailouts			
            G=tau*(Y-ZA*AB);
			DWLI=zetaDWL*delta*(1-ZR)*(ZM*qM*MI+ZA*qA*AI);
			bailout=F_eps_plus - (1-F_eps)*(WI - DWLI);
            BG=(1-tauL)*bailout/(qf-(1-tauL));
			Tvec=[shareB, shareI, shareS]' * (tauL*(BG+bailout));

            % saver consumption, multiplier
            Kbar=KBbar/shareHB;
            HS=Kbar*shareHS;
            BS = BI+BG;
            cS = (1-tau)*shareS*Y - qf*BS + BS/pi - nuK_budget.*p*HS - Tvec(3);
            Apay = zeta*ZA*(1+delta*(1-ZR)*qA);
            Mpay = zeta*(X + ZM*(1-delta+delta*ZR) + zeta*delta*ZM*(1-ZR)*qM);
            muS = max(0,qM + rbar*qA - betaS_g*(Mpay + rbar*Apay));
            
            % first-order conditions for borrowers
            % for p
            fx(1)= p - betaB_g*(rent + Del_BC_K)/( 1 - phiK*lamLTV - betaB_g*((1.0-ZK)*etaB + ZK*(1-nuK - (1-ZR)*phiK*lamLTV)) );
            % for kappabar
            fx(2)= kappabar - (1-OmM-rbar*OmA)*(1-delta*ZM*MB/(ZN*Mstar)) ...
                            + OmA*(rstar-rbar) ...
                            + p*phiK*lamLTV*(ZN*Kstar-ZK*KB)/(ZN*Mstar);
            % budget
            fx(3)= cB - ( (1-tau)*shareB*Y -ZA*(1-tau)*AB -ZM*(1-delta)*MB - delta*ZR*ZM*MB...
                           - p*SREO*KREO - rent*KREO - nuK_budget.*p*ZK*KB + ZR*ZN*Mstar - Tvec(1));
			% intermediary VF
			fx(4)= VI - (WI + CscrI); 	
            fx(5) = QK - (1-nuK - (1-ZR)*phiK*lamLTV)*p*KB;
            fx(6) = QA - ( (1-tau)*AB + delta*(1-ZR)*(OmA*AB) );
            fx(7) = QM - ( (1-delta + delta*ZR)*MB + delta*(1-ZR)*(OmM*MB) ); % Don't need Del_BC_M here, effects already included in OmM (DG)
            
 			% saver
            fx(8)= cy - xiL*cS/((1-xiL-xi)*BS);
            
            fx(9) = KB - KBbar/(1 + (1-ZK)/SREO);
            fx(10) = MB- (Lstar/(mu_G*pi-delta*ZM*(1-ZR)));
            fx(11) = AB - rstar*MB;
                                                            
            % intermediary and saver consumption
            HI=Kbar*shareHI;
            cI = F_eps*WI - F_eps_minus + qf*BI + (1-tau)*shareI*Y -delta*(1-ZR)*(ZM*qM*MI+ZA*qA*AI) ...
                    + p*(SREO-nuREO)*KREO + rent*KREO - pREO*(1-ZK)*KB - Lstar - nuK_budget.*p*HI - Tvec(2);            
            
            % check that goods market adds up
            Ycheck = cB + cI + cS + G + nuK_budget.*p*(ZK*KB + HI + HS) + nuREO*p*(KREO + (1-ZK)*KB) + (1-F_eps)*DWLI;            
            
            if print
                % print steady state values
                disp(' ');
                disp('Analytic steady state');
                disp('--- Rates ---');
                disp(['rstar: ',num2str(rstar)]);
                disp(['rD: ',num2str(rD)]);
                disp('--- Housing ---');
                disp(['KBbar: ',num2str(KBbar)]);
                disp(['KB: ',num2str(KB)]);
                disp(['KREO: ',num2str(KREO)]);                
                disp(['Kstar: ',num2str(Kstar)]);                
                disp(['p: ',num2str(p)]);
                disp(['pREO: ',num2str(pREO)]);
                disp(['lamLTV: ',num2str(lamLTV)]);
                disp('--- Debt ---');
                disp(['AB: ',num2str(AB)]);
                disp(['MB: ',num2str(MB)]);
                disp(['debt value: ',num2str(qM*MB+qA*AB)]);
                disp(['Mstar: ',num2str(Mstar)]);
                disp(['Lstar: ',num2str(Lstar)]);
                disp(['qM: ',num2str(qM)]);
                disp(['qA: ',num2str(qA)]);
                disp(['lamI: ',num2str(lamI)]);
                disp(['BI: ',num2str(BI)]);
                disp(['debt/GDP: ',num2str(qA*AB+qM*MB)]);
                disp(['debt/YB: ',num2str((qA*AB+qM*MB)/shareB)]);
                disp(['debt/HB: ',num2str(MB/(shareHB*p))]);
                disp(['VTI B: ',num2str(shareHB*p/shareB)]);
                disp(['debt/housing: ',num2str((qA*AB+qM*MB)/p)]);
                disp(['H exp: ',num2str(rent*KBbar/(shareB+rent*KB))]);
                disp('--- Liq Defs ---');
                disp(['PV(Del_BC_K): ',num2str( Del_BC_K/( 1 - phiK*lamLTV - betaB_g*ZK*(1-nuK - (1-ZR)*phiK*lamLTV) ) )]);
                disp(['PV(Del_BC_M): ',num2str( Del_BC_M/( 1 - phiI*lamI - betaB_g*zeta*delta*ZM*(1-ZR)) )]);
                disp(['PV(Del_BC_A): ',num2str( Del_BC_A/( 1 - phiI*lamI - betaB_g*zeta*delta*ZA*(1-ZR)) )]);
               disp('--- Borrower Risk ---');
                disp(['kappbar: ',num2str(kappabar)]);
                disp(['ZR: ',num2str(ZR)]);
                disp(['ZN: ',num2str(ZN)]);
                disp(['ZK: ',num2str(ZK)]);
                disp(['ZM: ',num2str(ZM)]);
                disp(['ZA: ',num2str(ZA)]);
                disp(['Drate ZN: ',num2str(1-ZN)]);
                disp(['X: ',num2str(X)]);
                disp('--- Intermediary Risk ---');
                disp(['1-F_eps: ',num2str(1-F_eps)]);
                disp(['DWL: ',num2str((1-F_eps)*DWLI)]);
                disp(['CscrI: ',num2str(CscrI)]);
                disp(['BG: ',num2str(BG)]);
                disp('--- Saver ---');
                disp(['cy: ',num2str(cy)]);
                disp(['muS: ',num2str(muS)]);                
                disp('--- Consumption ---');
                disp(['CB: ',num2str(cB)]);
                disp(['CI: ',num2str(cI)]);
                disp(['CS: ',num2str(cS)]);
                disp(['WI: ',num2str(WI)]);
                disp(['G: ',num2str(G)]);                
                disp(['maint: ',num2str(Y-cB-cI-cS-G)]);                
                disp(['Ycheck: ',num2str(Ycheck)]);                
            end
                      
            Sol=struct('qA',qA,...
                'qM',qM,...
                'qf',qf,...
                'p',p,...
                'pREO',pREO,...
                'cB',cB,...
                'cI',cI,...
                'cS',cS,...
                'ZR',ZR,... 
                'MS_g',0.1,...
                'lamLTV',lamLTV,...
                'lamI',lamI,...
                'muS',muS^(1/3),...
                'lamS',-0.3);

                
            rho=1-1/psi; 
            if psi==1 || 1-betaI*exp(rho*g)<0       
                VB=cB * exp(betaB * g / (1 - betaB) );
                VIHH=cI * exp(betaI * g / (1 - betaI) );
                VS=cS * exp(betaS * g / (1 - betaS) );
            else
                VB=cB * ((1-betaB)/(1-betaB*exp(rho*g)))^(1/rho);
                VIHH=cI * ((1-betaI)/(1-betaI*exp(rho*g)))^(1/rho);
                VS=cS * ((1-betaS)/(1-betaS*exp(rho*g)))^(1/rho);
            end
            
            frac_liq = theta * (1.0 - Z_vals(5)) / (1.0 - ZN);
            
            V=struct('cB',cB,...
                     'cI',cI,...
                     'cS',cS,...
                     'VB',VB,...
                     'VIHH',VIHH,...
                     'VS',VS,...
                     'OmA_B',OmA,...
                     'OmM_B',OmM,...
                     'OmK_B',OmK,...
                     'qA',qA,...
                     'qM',qM,...
                     'ZR',ZR,...
                     'p',p,...
                     'pREO',pREO,...
					 'CscrI',CscrI,...
                     'lamLTV',lamLTV,...
                     'BSnext',BS);
 
 
            nn = 1;
            
            ZN_str = Z_vals(nn); nn = nn + 1;
            ZK_str = Z_vals(nn); nn = nn + 1;
            ZM_str = Z_vals(nn); nn = nn + 1;
            ZA_str = Z_vals(nn); nn = nn + 1;
            
            ZN_liq = Z_vals(nn); nn = nn + 1;
            ZK_liq = Z_vals(nn); nn = nn + 1;
            ZM_liq = Z_vals(nn); nn = nn + 1;
            ZA_liq = Z_vals(nn); nn = nn + 1;                 
                 
            Add=struct('ZK',ZK,...
                         'ZN',ZN,...
                         'ZA',ZA,...
                         'ZM',ZM,...
                         'ZK_str',ZK_str,...
                         'ZN_str',ZN_str,...
                         'ZA_str',ZA_str,...
                         'ZM_str',ZM_str,...
                         'ZK_liq',ZK_liq,...
                         'ZN_liq',ZN_liq,...
                         'ZA_liq',ZA_liq,...
                         'ZM_liq',ZM_liq,...
                         'F_eps',F_eps,...
                         'F_eps_minus',F_eps_minus,...
                         'F_eps_plus',F_eps_plus,...
					     'VI', VI, ...		
                         'Mstar',Mstar,...
                         'Lstar',Lstar,...
                         'Kstar',Kstar,...
                         'AB_g',AB/mu_G,...
                         'MB_g',MB/mu_G,...
                         'BSnext',BS/G,...
                         'BInext',BI/G,...
                         'BGnext',BG/G,...
                         'cy',cy,...
                         'bailout',bailout,...
                         'X',X,...
                         'KREO_g',KREO/mu_G,...
                         'KB',KB,...
                         'MB',MB,...
                         'AB',AB,...
                         'OmA',OmA,...
                         'OmM',OmM,...
                         'sig2_om', sig2_om,...
                         'sig2_om_next',sig2_om,...
                         'rent',rent,...                        
                         'G',G);
            
                              
            State=struct('AB',AB,...
                         'MB',MB,...   
                         'KREO',KREO,...
                         'WI',WI,...
                         'BG',BG,...
                         'Y',Y,...
                         'G',G);
                       
            statsout=struct('KBbar',KBbar,...
                            'HI',HI,...
                            'HS',HS, ...
                            'p',p, ...
                            'ZN', ZN, ...
                            'frac_liq',frac_liq,...
                            'strat_pen',etaB*p*KBbar/VB);
            
            stvals=struct('Sol',Sol,...
                'V',V,...
                'Add',Add,...
                'State',State,...
                'statsout',statsout);  
            
        end
                      
        
        function [solguessvec,Vguessvec,V_names]=assignGuess(stv)

           solguess=struct('qA',log(stv.Sol.qA),...
                'qM',log(stv.Sol.qM),...
                'qf',log(stv.Sol.qf),...
                'p',log(stv.Sol.p),...
                'pREO',log(stv.Sol.pREO),...
                'cB',log(stv.Sol.cB),...
                'cI',log(stv.Sol.cI),...
                'cS',log(stv.Sol.cS),...
                'ZR',log(stv.Sol.ZR),...
                'MS_g',stv.Sol.MS_g,...
                'lamLTV',stv.Sol.lamLTV^(1/3),...
                'lamI',0.5,...
                'muS',stv.Sol.muS,...
                'lamS',stv.Sol.lamS);
            
            solguessvec=model.DSGEModel.structToVec(solguess);
                     
            Vguess=struct('cB',stv.V.cB,...
                'cI',stv.V.cI,...
                'cS',stv.V.cS,...
                'VB',stv.V.VB,...
                'VIHH',stv.V.VIHH,...
                'VS',stv.V.VS,...
                'exp_OmA_B',stv.V.OmA_B,...
                'exp_OmM_B',stv.V.OmM_B,...
                'exp_OmK_B',stv.V.OmK_B,...
                'qA',stv.V.qA,...
                'qM',stv.V.qM,...
                'ZR',stv.V.ZR,...
                'p',stv.V.p,...
                'pREO',stv.V.pREO,...
				'CscrI',stv.V.CscrI,...
                'lamLTV',stv.Sol.lamLTV^(1/3),...
                'BSnext',stv.V.BSnext);
            
            Vguessvec=model.DSGEModel.structToVec(Vguess);
            V_names=fieldnames(Vguess);
            
        end
                
        
        function [x,fx,exit,i]=tryOtherGuesses(fhand,gvec,options)
            % list of guesses                       

            gindex={[11,12],[11,12],[12,13],[10,12,13]};
            gvals={[-0.5,-0.5],[0.25,0.5],[-0.3,-0.5],[0.1,0.25,-0.25],};            
            
            for i=1:length(gindex)
                newguess=gvec;
                newguess(gindex{i})=gvals{i};
                [x,fx,exit]=fsolve(fhand,newguess,options);
%                if exit>=0
                if exit>0
                    break;
                end                                
            end
        end
                     
        
    end % of static methods
    
    
end

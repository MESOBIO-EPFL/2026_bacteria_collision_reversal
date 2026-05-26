%% Run a simulation with fixed and moving boundaries

function [bctAll,dtAll,fxBnd,mvBnd,dtBndAll,mp]=bct_simDynamicsFxMvBnd...
    (bctLn,nFa,vf,rto,prBndX,prBndY,...
    tsOrn,tsRev,tsCon,tsMmt,tsNmt,tsBnd,...
    clHrd,atMmt,mcAng,nmFrt,bndLn,...
    tsDiv,dvOn,...
    ttTm,tmSt,svPt)

%% Rnadomize number generator
rng('shuffle');

%% Set up model parameter structure array
addpath('./initial_configuration/')

mp=bct_parVarGenerator(bctLn,nFa,vf,rto,prBndX,prBndY,...
    tsOrn,tsRev,tsCon,tsMmt,tsNmt,tsBnd,...
    clHrd,atMmt,mcAng,nmFrt,bndLn,...
    tsDiv,dvOn,...
    ttTm,tmSt,svPt);

%% Creat initial configuration
[bctAll,dtAll,mp]=bct_iniConf(mp);
dtBndAll=dtAll;
fxBnd=bct_iniFxBnd(mp,2);
nePr=bct_neiPairRectangleAll(bctAll,mp);
nePrFxBnd=bct_neiPairFxBndAll(bctAll,fxBnd,mp);

rmpath('./initial_configuration/')

%% Initial relaxation of the system: 10\tau_R
[tmStOrg,tsNmtOrg,mtCellIdOrg]=deal(mp.tmSt,mp.tsNmt,mp.mtCellId);
mp.tmSt=0.001;
mp.tsNmt=1;
mp.mtCellId=ones(mp.nFa,1);

tmx=floor(10/mp.tmSt);

for tmc=1:tmx
    bctAll=bct_iterationFxBnd(bctAll,fxBnd,nePr,nePrFxBnd,mp);
    mxDis=bct_maximumDisplacement(bctAll,mp);
    if mxDis>0.48*mp.rngBf
        nePr=bct_neiPairRectangleAll(bctAll,mp);
        nePrFxBnd=bct_neiPairFxBndAll(bctAll,fxBnd,mp);
        mp.sdRngRef=bctAll(:,10:13);
    end
end

%% Reset model parameter values after initial relaxation
mp.tmSt=tmStOrg;
mp.tsNmt=tsNmtOrg;
mp.mtCellId=mtCellIdOrg;
bctAll(:,7:8)=0;

%% Dynamics simulation part

addpath('./initial_configuration/')
fxBnd=bct_iniFxBnd(mp,1);
mvBnd=bct_iniMvBnd(mp);
nePr=bct_neiPairRectangleAll(bctAll,mp);
nePrFxBnd=bct_neiPairFxBndAll(bctAll,fxBnd,mp);
nePrMvBnd=bct_neiPairMvBndAll(bctAll,mvBnd,mp);
rmpath('./initial_configuration/')

while floor(mp.curTm/mp.tmSt)<floor(mp.ttTm/mp.tmSt)
    [bctAll,mvBnd]=bct_iterationFxMvBnd(bctAll,fxBnd,mvBnd,...
        nePr,nePrFxBnd,nePrMvBnd,mp);

    bndEdLn=bct_boundaryEdgeLength(mvBnd,mp);
    if max(bndEdLn)>1.5*mp.bctLn || min(bndEdLn)<0.5*mp.bctLn
        mvBnd=bct_movingBoundaryRedistribute(mvBnd,bndEdLn,mp);
        nePrMvBnd=bct_neiPairMvBndAll(bctAll,mvBnd,mp);
    end
    mxDis=bct_maximumDisplacement(bctAll,mp);
    if mxDis>0.48*mp.rngBf
        nePr=bct_neiPairRectangleAll(bctAll,mp);
        nePrFxBnd=bct_neiPairFxBndAll(bctAll,fxBnd,mp);
        nePrMvBnd=bct_neiPairMvBndAll(bctAll,mvBnd,mp);
        mp.sdRngRef=bctAll(:,10:13);
    end

    if mp.dvOn==1
        bctAll(:,9)=bctAll(:,9)+(mp.bctLn+1)./mp.dvTm*mp.tmSt;
    end

    if mod(mp.itc,mp.dvChkTm)==1 && mp.dvOn==1
        [bctAll,mp,dvOcc]=bct_cellDivision(bctAll,mp);
        if dvOcc==1
            nePr=bct_neiPairRectangleAll(bctAll,mp);
            nePrFxBnd=bct_neiPairFxBndAll(bctAll,fxBnd,mp);
            nePrMvBnd=bct_neiPairMvBndAll(bctAll,mvBnd,mp);
            mp.sdRngRef=bctAll(:,10:13);
        end
    end

    mp.curTm=mp.curTm+mp.tmSt;
    mp.itc=mp.itc+1;

    if mp.svPt==1
        dtAll{mp.dtc}=bctAll;
        dtBndAll{mp.dtc}=mvBnd;
        mp.dtc=mp.dtc+1;
    elseif mod(mp.itc,mp.svPt)==1 && mp.svPt>1
        dtAll{mp.dtc}=bctAll;
        dtBndAll{mp.dtc}=mvBnd;
        mp.dtc=mp.dtc+1;
        mp.dtc
    end   
end

end
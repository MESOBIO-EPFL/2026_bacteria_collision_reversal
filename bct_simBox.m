%% Maze simulation function

function [dtAll,mp,bctAll,fxBnd]=bct_simBox...
    (bctLn,nFa,vf,rto,prBndX,prBndY,...
    tsOrn,tsRev,tsCon,tsMmt,tsNmt,tsBnd,...
    clHrd,atMmt,mcAng,nmFrt,bndLn,...
    tsDiv,dvOn,...
    ttTm,tmSt,svPt)

%% Set up model parameter structure array
rng('shuffle');

addpath('./initial_configuration/')

mp=bct_parVarGenerator(bctLn,nFa,vf,rto,prBndX,prBndY,...
    tsOrn,tsRev,tsCon,tsMmt,tsNmt,tsBnd,...
    clHrd,atMmt,mcAng,nmFrt,bndLn,...
    tsDiv,dvOn,...
    ttTm,tmSt,svPt);
mp.mzAng=pi/6;
mp.mzLn=sqrt(3*mp.nFa/mp.vf*(mp.bctLn+1)/sin(2*mp.mzAng)...
    +(25/2*cos(mp.mzAng)/sin(2*mp.mzAng))^2)-...
    (25/2*cos(mp.mzAng)/sin(2*mp.mzAng));

%% Creat initial configuration
[bctAll,dtAll,mp]=bct_iniConfMaze(mp);
fxBnd=bct_iniBoxFxBnd(mp,4);
nePr=bct_neiPairRectangleAll(bctAll,mp);
nePrFxBnd=bct_neiPairFxBndAll(bctAll,fxBnd,mp);

rmpath('./initial_configuration/')


%% Initial relaxation of the system: 10\tau_R
[tmStOrg,tsNmtOrg,mtCellIdOrg]=deal(mp.tmSt,mp.tsNmt,mp.mtCellId);
mp.tmSt=0.001;
mp.tsNmt=1;
mp.mtCellId=ones(mp.nFa,1); 
tmx=floor(5/mp.tmSt);

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

while floor(mp.curTm/mp.tmSt)<floor(mp.ttTm/mp.tmSt)
    bctAll=bct_iterationFxBnd(bctAll,fxBnd,nePr,nePrFxBnd,mp);
    mxDis=bct_maximumDisplacement(bctAll,mp);
    if mxDis>0.48*mp.rngBf
        nePr=bct_neiPairRectangleAll(bctAll,mp);
        nePrFxBnd=bct_neiPairFxBndAll(bctAll,fxBnd,mp);
        mp.sdRngRef=bctAll(:,10:13);
    end

    if mp.dvOn==1
        bctAll(:,9)=bctAll(:,9)+(mp.bctLn+1)./mp.dvTm*mp.tmSt;
    end

    if mod(mp.itc,mp.dvChkTm)==1 && mp.dvOn==1
        [bctAll,mp,dvOcc]=bct_cellDivision(bctAll,mp);
        if dvOcc==1
            nePr=bct_neiPairRectangleAll(bctAll,mp);
            mp.sdRngRef=bctAll(:,10:13);
        end
    end

    mp.curTm=mp.curTm+mp.tmSt;
    mp.itc=mp.itc+1;

    if mp.svPt==1
        dtAll{mp.dtc}=bctAll;
        mp.dtc=mp.dtc+1;
    elseif mod(mp.itc,mp.svPt)==1 && mp.svPt>1
        dtAll{mp.dtc}=bctAll;
        mp.dtc=mp.dtc+1;
        mp.dtc
    end   
end

end
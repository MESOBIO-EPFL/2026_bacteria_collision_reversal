%% Initialize workspace
clearvars;
close all;

%% Input parameter values
[bctLn,nFa,vf,rto,prBndX,prBndY]=deal(1.5,100,0.7,2,0,1);
[tsOrn,tsRev,tsCon,tsMmt,tsNmt,tsBnd]=deal(1,-1,1,10,1,3);
[clHrd,atMmt,mcAng,nmFrt,bndLn]=deal(100,0.2,pi/6,0,10);
[tsDiv,dvOn]=deal(10,0);
[ttTm,tmSt,svPt]=deal(100,0.001,500);

%% Run simulations
[bctAll,dtAll,fxBnd,mvBnd,dtBndAll,mp]=bct_simDynamicsFxMvBnd...
    (bctLn,nFa,vf,rto,prBndX,prBndY,...
    tsOrn,tsRev,tsCon,tsMmt,tsNmt,tsBnd,...
    clHrd,atMmt,mcAng,nmFrt,bndLn,...
    tsDiv,dvOn,...
    ttTm,tmSt,svPt);

%% Set up file name
csNm=sprintf(['bctFxBnd_n%d_vf%s_bctLn%s_atMmt%s',...
    '_tsCon%s_nmFrt%s_dvOn%d'],...
    mp.nFa,num2str(mp.vf),num2str(mp.bctLn),...
    num2str(mp.atMmt),num2str(mp.tsCon),...
    num2str(mp.nmFrt),mp.dvOn);
csNm=strrep(csNm,'.','d');
odr='./data/';
rpc=1;

%% Save the data file
bct_saveDataFxMvBnd(dtAll,dtBndAll,fxBnd,mp,csNm,odr,rpc);

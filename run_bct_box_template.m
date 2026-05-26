%% Box simulation

%% Initialize workspace
clearvars;
close all;

%% Input parameter values
[bctLn,nFa,vf,rto,prBndX,prBndY]=deal(1.5,200,0.4,1,0,0);
[tsOrn,tsRev,tsCon,tsMmt,tsNmt,tsBnd]=deal(1,-1,0.1,10,1,10);
[clHrd,atMmt,mcAng,nmFrt,bndLn]=deal(50,0.2,pi/6,0,10);
[tsDiv,dvOn]=deal(10,0);
[ttTm,tmSt,svPt]=deal(50,0.001,1000);

%% Run simulation
[dtAll,mp,bctAll,fxBnd]=bct_simBox...
    (bctLn,nFa,vf,rto,prBndX,prBndY,...
    tsOrn,tsRev,tsCon,tsMmt,tsNmt,tsBnd,...
    clHrd,atMmt,mcAng,nmFrt,bndLn,...
    tsDiv,dvOn,...
    ttTm,tmSt,svPt);

%% Set up file name
csNm=sprintf(['bctBox_n%d_vf%s_bctLn%s_atMmt%s',...
    '_tsCon%s_nmFrt%s_dvOn%d'],...
    mp.nFa,num2str(mp.vf),num2str(mp.bctLn),...
    num2str(mp.atMmt),num2str(mp.tsCon),...
    num2str(mp.nmFrt),mp.dvOn);
csNm=strrep(csNm,'.','d');

%% Save the data file
rpc=1;
odr='./data/';
bct_saveDataFxBnd(dtAll,fxBnd,mp,csNm,odr,rpc);


%% Initialize workspace
clearvars;
close all;

%% Input parameter values
[bctLn,nFa,vf,rto,prBndX,prBndY]=deal(2,50,0.7,1,1,1);
[tsOrn,tsRev,tsCon,tsMmt,tsNmt,tsBnd]=deal(1,-1,0,100,1,10);
[clHrd,atMmt,mcAng,nmFrt,bndLn]=deal(50,0.1,pi/6,0,10);
[tsDiv,dvOn]=deal(10,0);
[ttTm,tmSt,svPt]=deal(50,0.0002,2500);

%% Run a simulation
[bctAll,dtAll,mp]=bct_simDynamics...
    (bctLn,nFa,vf,rto,prBndX,prBndY,...
    tsOrn,tsRev,tsCon,tsMmt,tsNmt,tsBnd,...
    clHrd,atMmt,mcAng,nmFrt,bndLn,...
    tsDiv,dvOn,...
    ttTm,tmSt,svPt);

%% Set up file name
rpc=1;
csNm=sprintf(['bctDyn_n%d_vf%s_bctLn%s_atMmt%s',...
    '_tsCon%s_nmFrt%s_dvOn%d'],...
    mp.nFa,num2str(mp.vf),num2str(mp.bctLn),...
    num2str(mp.atMmt),num2str(mp.tsCon),...
    num2str(mp.nmFrt),mp.dvOn);
csNm=strrep(csNm,'.','d');
odr='./data/';

%% Save the data file
bct_saveData(dtAll,odr,mp,csNm,rpc);


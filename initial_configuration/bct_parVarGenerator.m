%% Create a structure variable, store all relevant parameters

function mp=bct_parVarGenerator...
    (bctLn,nFa,vf,rto,prBndX,prBndY,...
    tsOrn,tsRev,tsCon,tsMmt,tsNmt,tsBnd,...
    clHrd,atMmt,mcAng,nmFrt,bndLn,...
    tsDiv,dvOn,...
    ttTm,tmSt,svPt)

%% Set up geometric parameters
mp=struct;
[mp.bctLn,mp.nFa,mp.vf]=deal(bctLn,nFa,vf);
[mp.rto,mp.prBndX,mp.prBndY]=deal(rto,prBndX,prBndY);
mp.lx=sqrt(nFa*(3*bctLn/2+1/2+pi/4)/rto/vf);
mp.ly=rto*mp.lx;

%% Set up model parameters
[mp.tsOrn,mp.tsRev,mp.tsCon,mp.tsMmt]=deal(tsOrn,tsRev,tsCon,tsMmt);
[mp.tsNmt,mp.tsBnd]=deal(tsNmt,tsBnd);

[mp.clHrd,mp.atMmt,mp.mcAng]=deal(clHrd,atMmt,mcAng);
[mp.nmFrt,mp.bndLn]=deal(nmFrt,bndLn);

nmCellMx=floor(mp.nmFrt*mp.nFa);
mp.mtCellId=ones(mp.nFa,1);
mp.tmSclFrt=ones(mp.nFa,1);

if nmCellMx>0
    nmCellId=randperm(mp.nFa,nmCellMx);
    mp.mtCellId(nmCellId)=0;
    mp.tmSclFrt(nmCellId)=tsNmt;
end

[mp.ttTm,mp.tmSt,mp.svPt]=deal(ttTm,tmSt,svPt);
[mp.itc,mp.curTm,mp.dtc]=deal(0,0,1);

mp.rngBf=0.3;

[mp.tsDiv,mp.dvOn]=deal(tsDiv,dvOn);
mp.vfMx=0.6;
mp.dvChkTm=5;

end
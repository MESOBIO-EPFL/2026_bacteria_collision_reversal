%% Individual iteration step (Periodic)

function bctAllN=bct_iteration(bctAll,nePr,mp)

%% Compute forces (passive and active) and moments (passive and active)
[frc,mnt,bctAllN]=bct_forceComputeAll(bctAll,nePr,mp);
frcAct=bctAllN(:,4).*mp.mtCellId.*[cos(bctAllN(:,3)),sin(bctAllN(:,3))];

if mp.atMmt>0
    bctAllN(:,5)=bctAllN(:,5).*(1-mp.tmSt/mp.tsMmt)+mp.atMmt...
        .*mp.mtCellId.*sqrt(2/mp.tsMmt).*sqrt(mp.tmSt)...
        .*normrnd(0,1,mp.nFa,1);
end

%% Update cell position and orientation
bctAllN(:,1:2)=bctAllN(:,1:2)+1./mp.tmSclFrt.*(frc+frcAct).*mp.tmSt;
bctAllN(:,3)=bctAllN(:,3)+1/mp.tsOrn./mp.tmSclFrt.*...
    (mnt+bctAllN(:,5)).*mp.tmSt;

if mp.prBndX==1
    bctAllN(:,1)=mod(bctAllN(:,1),mp.lx);
end

if mp.prBndY==1
    bctAllN(:,2)=mod(bctAllN(:,2),mp.ly);
end

%% Apply spontaneous reversal for each cell
if any(bctAllN(:,6)>mp.tsRev) && mp.tsRev>0
    bctAllN((bctAllN(:,6)>mp.tsRev),4)=...
        (-1)*bctAllN((bctAllN(:,6)>mp.tsRev),4);
    bctAllN(bctAllN(:,6)>mp.tsRev,7)=bctAllN(bctAllN(:,6)>mp.tsRev,7)+1;
    bctAllN(bctAllN(:,6)>mp.tsRev,6)=0;
end

%% Update cell range
bctAllN(:,10:13)=bct_cellRangeCompute...
    (bctAllN(:,1:2),bctAllN(:,3),bctAllN(:,9),mp.rngBf);

end
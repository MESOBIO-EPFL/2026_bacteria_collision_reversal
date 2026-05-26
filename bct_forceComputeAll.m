%% Compute force for each cell

function [frc,mnt,bctAllN]=bct_forceComputeAll(bctAll,nePr,mp)

%% Define output variables
bctAllN=bctAll;
frc=zeros(mp.nFa,2);
mnt=zeros(mp.nFa,1);

%% Find all crosing cells, separate into crossing/non-crossing groups
[cross,parVal]=bct_lineCrossAll(bctAllN,nePr,mp);

%% Seperate crossing pairs(nePrCr) and non-crossing pairs(nePrNcr)
nePrCr=nePr(cross==1,:);
pvCr=parVal(cross==1,:);
nePrNcr=nePr(cross==0,:);

if ~isempty(nePrCr)
    for prc=1:size(nePrCr,1)
        %% first, locally correct cell positions
        lcCen=bct_crdLocal(bctAllN(nePrCr(prc,:),1:2),...
            mp.lx,mp.ly,mp.prBndX,mp.prBndY);
        lcOrn=bctAllN(nePrCr(prc,:),3);
        lcLen=bctAllN(nePrCr(prc,:),9);

        %% If there is crossing, compute force
        cenVec=lcCen(1,:)-lcCen(2,:);
        if abs(pvCr(prc,1)-0.5)<abs(pvCr(prc,2)-0.5)
            frVec=[-sin(lcOrn(1)),cos(lcOrn(1))];
        else
            frVec=[-sin(lcOrn(2)),cos(lcOrn(2))];
        end

        if dot(cenVec,frVec)<0
            frVec=-frVec;
        end
        frVec=frVec*mp.clHrd;

        frc(nePrCr(prc,1),:)=frc(nePrCr(prc,1),:)+frVec;
        frc(nePrCr(prc,2),:)=frc(nePrCr(prc,2),:)-frVec;

        %% Update moment
        armVec=[lcLen(1)*(pvCr(prc,1)-0.5)*cos(lcOrn(1)),...
            lcLen(1)*(pvCr(prc,1)-0.5)*sin(lcOrn(1));...
            lcLen(2)*(pvCr(prc,2)-0.5)*cos(lcOrn(2)),...
            lcLen(2)*(pvCr(prc,2)-0.5)*sin(lcOrn(2))];

        mnt(nePrCr(prc,1))=mnt(nePrCr(prc,1))+...
            armVec(1,1)*frVec(2)-armVec(1,2)*frVec(1);
        mnt(nePrCr(prc,2))=mnt(nePrCr(prc,2))-...
            armVec(2,1)*frVec(2)+armVec(2,2)*frVec(1);    

        %% Contact reversal calculation
        if mp.tsCon>0
            angDff=mod(lcOrn(2)-lcOrn(1),pi);
            if angDff>pi/2
                angDff=pi-angDff;
            end
            if angDff>mp.mcAng
                cnPrb=(1-0.5*((angDff-mp.mcAng)/(pi/2-mp.mcAng)))...
                    ^(mp.tmSt/mp.tsCon);
                if rand()>cnPrb
                    if abs(pvCr(prc,1)-0.5)>abs(pvCr(prc,2)-0.5)
                        bctAllN(nePrCr(prc,1),4)=(-1)*...
                            bctAllN(nePrCr(prc,1),4);
                        bctAllN(nePrCr(prc,1),6)=0;
                        bctAllN(nePrCr(prc,1),8)=...
                            bctAllN(nePrCr(prc,1),8)+1;
                    else
                        bctAllN(nePrCr(prc,2),4)=(-1)*...
                            bctAllN(nePrCr(prc,2),4);
                        bctAllN(nePrCr(prc,2),6)=0;
                        bctAllN(nePrCr(prc,2),8)=...
                            bctAllN(nePrCr(prc,2),8)+1;
                    end
                end
            end
        end
    end
end

if ~isempty(nePrNcr)
    [minDis,minPar,minCrd]=bct_minimumDistanceAll(bctAllN,nePrNcr,mp);

    valChk=(minDis>0 & minDis<=1);
    minDisVal=minDis(valChk);
    minParVal=minPar(valChk,:);
    minCrdVal=minCrd(valChk,:);
    nePrNcrVal=nePrNcr(valChk,:);

    if ~isempty(minCrdVal)

        frVec=minCrdVal(:,1:2)-minCrdVal(:,3:4);
        frVec=frVec./(sqrt(sum(frVec.^2,2)));
        frVec=mp.clHrd.*frVec.*(1-minDisVal);

        armVecP=[bctAllN(nePrNcrVal(:,1),9).*(minParVal(:,1)-0.5).*...
            cos(bctAllN(nePrNcrVal(:,1),3)),...
            bctAllN(nePrNcrVal(:,1),9).*(minParVal(:,1)-0.5).*...
            sin(bctAllN(nePrNcrVal(:,1),3))];
        armVecQ=[bctAllN(nePrNcrVal(:,2),9).*(minParVal(:,2)-0.5).*...
            cos(bctAllN(nePrNcrVal(:,2),3)),...
            bctAllN(nePrNcrVal(:,2),9).*(minParVal(:,2)-0.5).*...
            sin(bctAllN(nePrNcrVal(:,2),3))];

        for prc=1:numel(nePrNcrVal)/2
            frc(nePrNcrVal(prc,1),:)=frc(nePrNcrVal(prc,1),:)+frVec(prc,:);
            frc(nePrNcrVal(prc,2),:)=frc(nePrNcrVal(prc,2),:)-frVec(prc,:);

            mnt(nePrNcrVal(prc,1))=mnt(nePrNcrVal(prc,1))+...
                armVecP(prc,1).*frVec(prc,2)-armVecP(prc,2).*frVec(prc,1);
            mnt(nePrNcrVal(prc,2))=mnt(nePrNcrVal(prc,2))+...
                armVecQ(prc,1).*(-frVec(prc,2))+armVecQ(prc,2).*frVec(prc,1);
        end

        %% Compute contact reversal
        if mp.tsCon>0
            angDff=mod(bctAllN(nePrNcrVal(:,2),3)-...
                bctAllN(nePrNcrVal(:,1),3),pi);
            angDff(angDff>pi/2)=pi-angDff(angDff>pi/2);

            minParAng=minParVal(angDff>mp.mcAng,:);
            nePrNcrAng=nePrNcrVal(angDff>mp.mcAng,:);
            angDff=angDff(angDff>mp.mcAng);
            cnPrb=(1-0.5.*((angDff-mp.mcAng)./(pi/2-mp.mcAng)))...
                        .^(mp.tmSt/mp.tsCon);

            valFst=((abs(minParAng(:,1)-0.5))>(abs(minParAng(:,2)-0.5)));
            valSnd=((abs(minParAng(:,1)-0.5))<=(abs(minParAng(:,2)-0.5)));

            amx=numel(nePrNcrAng)/2;
            cnRevId=zeros(amx,1);
            cnRevId(valFst)=nePrNcrAng(valFst,1);
            cnRevId(valSnd)=nePrNcrAng(valSnd,2);

            rndNm=rand(amx,1);
            cnRevTrue=find(rndNm>cnPrb);

            if ~isempty(cnRevTrue)
                for prc=1:numel(cnRevTrue)
                    bctAllN(cnRevId(cnRevTrue(prc)),4)=(-1)*...
                        bctAllN(cnRevId(cnRevTrue(prc)),4);
                    bctAllN(cnRevId(cnRevTrue(prc)),6)=0;
                    bctAllN(cnRevId(cnRevTrue(prc)),8)=...
                        bctAllN(cnRevId(cnRevTrue(prc)),8)+1;
                end
            end
        end
    end
end

%% Update reversal time values
if mp.tsRev>0
    bctAllN(:,6)=bctAllN(:,6)+mp.tmSt;
end

end


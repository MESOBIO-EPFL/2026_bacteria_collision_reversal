%% Compute force for each cell

function [frcBnd,mntBnd,mvBndFrc,bctAllN]=bct_forceComputeMvBndAll...
    (bctAll,mvBnd,nePrMvBnd,mp)

%% Define output variables
bctAllN=bctAll;
frcBnd=zeros(mp.nFa,2);
mntBnd=zeros(mp.nFa,1);
mvBndFrc=zeros(size(mvBnd,1),2);
mvBndc=[mvBnd;mvBnd(1,:)];

if ~isempty(nePrMvBnd)
    %% Find all crosing cells, separate into crossing/non-crossing groups
    [cross,parVal]=bct_lineCrossMvBndAll(bctAllN,mvBnd,nePrMvBnd,mp);
    bndPtIn=mvBndc(1:end-1,:);
    bndPtFn=mvBndc(2:end,:);

    if mp.prBndX==1
        bndPtFn(:,1)=bndPtFn(:,1)-mp.lx.*...
            floor((bndPtFn(:,1)-bndPtIn(:,1))./mp.lx+1/2);
    end
    if mp.prBndY==1
        bndPtFn(:,2)=bndPtFn(:,2)-mp.ly.*...
            floor((bndPtFn(:,2)-bndPtIn(:,2))./mp.ly+1/2);
    end

    tnVec=bndPtFn-bndPtIn;
    tnVec=tnVec./sqrt(sum(tnVec.^2,2));
    nmVec=[-tnVec(:,2),tnVec(:,1)];

    %% Seperate crossing and non-crossing cell pairs
    nePrCr=nePrMvBnd(cross==1,:);
    pvCr=parVal(cross==1,:);
    nePrNcr=nePrMvBnd(cross==0,:);

    if ~isempty(nePrCr)
        for prc=1:size(nePrCr,1)
           %% first, locally correct cell positions
            lcOrn=bctAllN(nePrCr(prc,2),3);
            lcLen=bctAllN(nePrCr(prc,2),9);
            lcBnc=nePrCr(prc,1);

            %% If there is crossing, compute force
            frVec=nmVec(lcBnc,:)*mp.clHrd;
            frcBnd(nePrCr(prc,2),:)=frcBnd(nePrCr(prc,2),:)+frVec;

            %% Update moment
            armVec=[lcLen*(pvCr(prc,2)-0.5)*cos(lcOrn),...
                lcLen*(pvCr(prc,2)-0.5)*sin(lcOrn)];

            mntBnd(nePrCr(prc,2))=mntBnd(nePrCr(prc,2))+...
                armVec(1)*frVec(2)-armVec(2)*frVec(1);    

            %% Compute force to boundary vertices   
            mvBndFrc(lcBnc,:)=mvBndFrc(lcBnc,:)-(1-pvCr(prc,1))*frVec;
            if lcBnc==size(mvBnd,1)
                mvBndFrc(1,:)=mvBndFrc(1,:)-pvCr(prc,1)*frVec;
            else
                mvBndFrc(lcBnc+1,:)=mvBndFrc(lcBnc+1,:)-pvCr(prc,1)*frVec;
            end
        end
    end
    
    if ~isempty(nePrNcr)
        [minDis,minPar,minCrd]=bct_minimumDistanceMvBndAll...
            (bctAll,mvBnd,nePrNcr,mp);
    
        valChk=(minDis>0 & minDis<=0.5);
        minDisVal=minDis(valChk);
        minParVal=minPar(valChk,:);
        minCrdVal=minCrd(valChk,:);
        nePrNcrVal=nePrNcr(valChk,:);
    
        if ~isempty(minCrdVal)
   
            %% Compute force and moment to cells from boundary
            frVec=nmVec(nePrNcrVal(:,1),:);
            frVec=mp.clHrd.*frVec.*(1-2*minDisVal);    
            
            armVecQ=[bctAll(nePrNcrVal(:,2),9).*(minParVal(:,2)-0.5).*...
                cos(bctAll(nePrNcrVal(:,2),3)),...
                bctAll(nePrNcrVal(:,2),9).*(minParVal(:,2)-0.5).*...
                sin(bctAll(nePrNcrVal(:,2),3))];
    
    
            for prc=1:numel(nePrNcrVal)/2
                frcBnd(nePrNcrVal(prc,2),:)=...
                    frcBnd(nePrNcrVal(prc,2),:)+frVec(prc,:);
                mntBnd(nePrNcrVal(prc,2))=mntBnd(nePrNcrVal(prc,2))+...
                    armVecQ(prc,1).*(frVec(prc,2))...
                    -armVecQ(prc,2).*frVec(prc,1);
            end

            %% Compute force to boundary
            for prc=1:numel(nePrNcrVal)/2
                mvBndFrc(nePrNcrVal(prc,1),:)=...
                    mvBndFrc(nePrNcrVal(prc,1),:)-...
                    (1-minParVal(prc,1))*frVec(prc,:);
                if nePrNcrVal(prc,1)==size(mvBnd,1)
                    mvBndFrc(1,:)=mvBndFrc(1,:)-...
                        (minParVal(prc,1))*frVec(prc,:);
                else
                    mvBndFrc(nePrNcrVal(prc,1)+1,:)=...
                        mvBndFrc(nePrNcrVal(prc,1)+1,:)-...
                        (minParVal(prc,1))*frVec(prc,:);
                end
            end
        end
    end
end

end


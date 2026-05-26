%% Compute force for each cell

function [frcBnd,mntBnd,bctAllN]=bct_forceComputeFxBndAll...
    (bctAll,fxBnd,nePrFxBnd,mp)

%% Define output variables
bctAllN=bctAll;
frcBnd=zeros(mp.nFa,2);
mntBnd=zeros(mp.nFa,1);

if ~isempty(nePrFxBnd)
    %% Find all crosing cells, separate into crossing/non-crossing groups
    [cross,parVal]=bct_lineCrossFxBndAll(bctAllN,fxBnd,nePrFxBnd,mp);

    %% Seperate crossing and non-crossing cell pairs
    nePrCr=nePrFxBnd(cross==1,:);
    pvCr=parVal(cross==1,:);
    nePrNcr=nePrFxBnd(cross==0,:);

    if ~isempty(nePrCr)
        for prc=1:size(nePrCr,1)
           %% first, locally correct cell positions
            lcOrn=bctAllN(nePrCr(prc,2),3);
            lcLen=bctAllN(nePrCr(prc,2),9);

            %% If there is crossing, compute force
            frVec=fxBnd(nePrCr(prc,1),5:6)*mp.clHrd;
            frcBnd(nePrCr(prc,2),:)=frcBnd(nePrCr(prc,2),:)+frVec;

            %% Update moment
            armVec=[lcLen*(pvCr(prc,2)-0.5)*cos(lcOrn),...
                lcLen*(pvCr(prc,2)-0.5)*sin(lcOrn)];

            mntBnd(nePrCr(prc,2))=mntBnd(nePrCr(prc,2))+...
                armVec(1)*frVec(2)-armVec(2)*frVec(1);    
            %% Contact reversal calculation
            plVec=bctAllN(nePrCr(prc,2),4).*...
                [cos(bctAllN(nePrCr(prc,2),3)),...
                sin(bctAllN(nePrCr(prc,2),3))];
            if dot(fxBnd(nePrCr(prc,1),5:6),plVec)<0
                bctAllN(nePrCr(prc,2),4)=-bctAllN(nePrCr(prc,2),4);
            end
        end
    end
    
    if ~isempty(nePrNcr)
        [minDis,minPar,minCrd]=bct_minimumDistanceFxBndAll...
            (bctAll,fxBnd,nePrNcr,mp);
    
        valChk=(minDis>0 & minDis<=0.5);
        minDisVal=minDis(valChk);
        minParVal=minPar(valChk,:);
        minCrdVal=minCrd(valChk,:);
        nePrNcrVal=nePrNcr(valChk,:);
    
        if ~isempty(minCrdVal)
    
            frVec=fxBnd(nePrNcrVal(:,1),5:6);
            frVec=mp.clHrd.*frVec.*(1-2*minDisVal);
    
            
            armVecQ=[bctAll(nePrNcrVal(:,2),9).*(minParVal(:,2)-0.5).*...
                cos(bctAll(nePrNcrVal(:,2),3)),...
                bctAll(nePrNcrVal(:,2),9).*(minParVal(:,2)-0.5).*...
                sin(bctAll(nePrNcrVal(:,2),3))];
    
    
            for prc=1:numel(nePrNcrVal)/2
                frcBnd(nePrNcrVal(prc,2),:)=...
                    frcBnd(nePrNcrVal(prc,2),:)+frVec(prc,:);
                mntBnd(nePrNcrVal(prc,2))=mntBnd(nePrNcrVal(prc,2))+...
                    armVecQ(prc,1).*(frVec(prc,2))-armVecQ(prc,2).*frVec(prc,1);
            end

            %% Compute contact reversal
            plVec=bctAllN(nePrNcrVal(:,2),4).*...
                [cos(bctAllN(nePrNcrVal(:,2),3)),...
                sin(bctAllN(nePrNcrVal(:,2),3))];
            dtPrd=fxBnd(nePrNcrVal(:,1),5).*plVec(:,1)+...
                fxBnd(nePrNcrVal(:,1),6).*plVec(:,2);
            bctAllN(nePrNcrVal(dtPrd<0,2),4)=...
                -bctAllN(nePrNcrVal(dtPrd<0,2),4);
        end
    end
end

end


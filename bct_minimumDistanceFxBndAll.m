%% Compute minimum distance between two segments that do not cross

function [minDis,minPar,minCrd]=bct_minimumDistanceFxBndAll...
    (bctAll,bndAll,nePrNcr,mp)

nmx=numel(nePrNcr);

pEndIn=bndAll(nePrNcr(:,1),1:2);
pEndFn=bndAll(nePrNcr(:,1),3:4);
pVec=pEndFn-pEndIn;
pMagSq=sum(pVec.^2,2);


cnCrd=bctAll(nePrNcr(:,2),1:2);

qEndIn=[cnCrd(:,1)-bctAll(nePrNcr(:,2),9)/2.*cos(bctAll(nePrNcr(:,2),3)),...
    cnCrd(:,2)-bctAll(nePrNcr(:,2),9)/2.*sin(bctAll(nePrNcr(:,2),3))];
qEndFn=[cnCrd(:,1)+bctAll(nePrNcr(:,2),9)/2.*cos(bctAll(nePrNcr(:,2),3)),...
    cnCrd(:,2)+bctAll(nePrNcr(:,2),9)/2.*sin(bctAll(nePrNcr(:,2),3))];
qVec=qEndFn-qEndIn;
qMagSq=sum(qVec.^2,2);


pvNcr=[sum((pEndIn-qEndIn).*qVec,2)./qMagSq,...
    sum((pEndFn-qEndIn).*qVec,2)./qMagSq,...
    sum((qEndIn-pEndIn).*pVec,2)./pMagSq,...
    sum((qEndFn-pEndIn).*pVec,2)./pMagSq];

pvNcr(pvNcr<0)=0;
pvNcr(pvNcr>1)=1;

pvP=[repmat([0,1],nmx/2,1),pvNcr(:,3:4)];
pvQ=[pvNcr(:,1:2),repmat([0,1],nmx/2,1)];

dis=zeros(nmx/2,1);
for csc=1:4
    dis(:,csc)=sqrt(sum(((pEndIn+pvP(:,csc).*pVec)-...
        (qEndIn+pvQ(:,csc).*qVec)).^2,2));
end

[minDis,minIdx]=min(dis,[],2);
idx=sub2ind(size(pvP),(1:nmx/2).',minIdx);
minPar=[pvP(idx),pvQ(idx)];
minCrd=[pEndIn+minPar(:,1).*pVec,qEndIn+minPar(:,2).*qVec];

end
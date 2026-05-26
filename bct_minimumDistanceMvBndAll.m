%% Compute minimum distance between two segments that do not cross

function [minDis,minPar,minCrd]=bct_minimumDistanceMvBndAll...
    (bctAll,mvBnd,nePrNcr,mp)

nmx=numel(nePrNcr);
mvBndc=[mvBnd;mvBnd(1,:)];

pEndIn=mvBndc(nePrNcr(:,1),:);
pEndFn=mvBndc(nePrNcr(:,1)+1,:);

if mp.prBndX==1
    pEndFn(:,1)=pEndFn(:,1)-mp.lx.*...
        floor((pEndFn(:,1)-pEndIn(:,1))./mp.lx+1/2);
end
if mp.prBndY==1
    pEndFn(:,2)=pEndFn(:,2)-mp.ly.*...
        floor((pEndFn(:,2)-pEndIn(:,2))./mp.ly+1/2);
end

pVec=pEndFn-pEndIn;
pMagSq=sum(pVec.^2,2);


cnCrd=bctAll(nePrNcr(:,2),1:2);
if mp.prBndX==1
    cnCrd(:,1)=cnCrd(:,1)-mp.lx.*...
        floor((cnCrd(:,1)-pEndIn(:,1))./mp.lx+1/2);
end
if mp.prBndY==1
    cnCrd(:,2)=cnCrd(:,2)-mp.ly.*...
        floor((cnCrd(:,2)-pEndIn(:,2))./mp.ly+1/2);
end

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
%% Compute minimum distance between two segments that do not cross

function [minDis,minPar,minCrd]=bct_minimumDistanceAll(bctAll,nePrNcr,mp)

nePrNcrVec=(nePrNcr.');
nePrNcrVec=nePrNcrVec(:);
nmx=numel(nePrNcrVec);

%% locally adjust particle positions under periodic boundary condition
cnCrd=bctAll(nePrNcrVec,1:2);
if mp.prBndX==1
    cnCrd(2:2:nmx,1)=cnCrd(2:2:nmx,1)-mp.lx.*...
        floor((cnCrd(2:2:nmx,1)-cnCrd(1:2:nmx,1))./mp.lx+1/2);
end
if mp.prBndY==1
    cnCrd(2:2:nmx,2)=cnCrd(2:2:nmx,2)-mp.ly.*...
        floor((cnCrd(2:2:nmx,2)-cnCrd(1:2:nmx,2))./mp.ly+1/2);
end

%% Find coefficients for distance calculation
pEndIn=[cnCrd(1:2:nmx,1)-bctAll(nePrNcr(:,1),9)/2.*...
    cos(bctAll(nePrNcr(:,1),3)),...
    cnCrd(1:2:nmx,2)-bctAll(nePrNcr(:,1),9)/2.*...
    sin(bctAll(nePrNcr(:,1),3))];
pEndFn=[cnCrd(1:2:nmx,1)+bctAll(nePrNcr(:,1),9)/2.*...
    cos(bctAll(nePrNcr(:,1),3)),...
    cnCrd(1:2:nmx,2)+bctAll(nePrNcr(:,1),9)/2.*...
    sin(bctAll(nePrNcr(:,1),3))];
pVec=pEndFn-pEndIn;
pMagSq=sum(pVec.^2,2);

qEndIn=[cnCrd(2:2:nmx,1)-bctAll(nePrNcr(:,2),9)/2.*...
    cos(bctAll(nePrNcr(:,2),3)),...
    cnCrd(2:2:nmx,2)-bctAll(nePrNcr(:,2),9)/2.*...
    sin(bctAll(nePrNcr(:,2),3))];
qEndFn=[cnCrd(2:2:nmx,1)+bctAll(nePrNcr(:,2),9)/2.*...
    cos(bctAll(nePrNcr(:,2),3)),...
    cnCrd(2:2:nmx,2)+bctAll(nePrNcr(:,2),9)/2.*...
    sin(bctAll(nePrNcr(:,2),3))];
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
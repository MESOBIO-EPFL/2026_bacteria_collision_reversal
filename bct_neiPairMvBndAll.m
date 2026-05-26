%% Find all neighbor pairs, using disk radius

function nePrMvBnd=bct_neiPairMvBndAll(bctAll,mvBnd,mp)

%% Copy cells for finding all neighbors with fixed boundary edges
[~,~,~,sdRngCp,sdIdCp]=...
    bct_seedCopy(bctAll(:,1:2),bctAll(:,3),...
    bctAll(:,9),bctAll(:,10:13),...
    mp.lx,mp.ly,mp.prBndX,mp.prBndY,max(bctAll(:,9))+1);

nePrMvBnd=cell(size(mvBnd,1),1);
mvBndc=[mvBnd;mvBnd(1,:)];

%% Check box crossing for each fixed boundary
for edc=1:size(mvBnd,1)
    bndCrd=bct_crdLocal(mvBndc(edc:edc+1,:),...
        mp.lx,mp.ly,mp.prBndX,mp.prBndY);
    valCs=((sdRngCp(:,1)<max(bndCrd(:,1))) & ...
        (sdRngCp(:,2)>min(bndCrd(:,1))) & ...
        (sdRngCp(:,3)<max(bndCrd(:,2))) & ...
        (sdRngCp(:,4)>min(bndCrd(:,2))));
    neId=sdIdCp(valCs);
    nePrMvBnd{edc}=[ones(numel(neId),1)*edc,neId];
end

%% Create a pair of potential neighboring fixed boundary and cells
nePrMvBnd=nePrMvBnd(~cellfun('isempty',nePrMvBnd));
nePrMvBnd=cell2mat(nePrMvBnd);
nePrMvBnd=unique(nePrMvBnd,'rows');

end
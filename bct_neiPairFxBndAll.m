%% Find all neighbor pairs, using disk radius

function nePrFxBnd=bct_neiPairFxBndAll(bctAll,fxBnd,mp)

%% Copy cells for finding all neighbors with fixed boundary edges
[~,~,~,sdRngCp,sdIdCp]=...
    bct_seedCopy(bctAll(:,1:2),bctAll(:,3),...
    bctAll(:,9),bctAll(:,10:13),...
    mp.lx,mp.ly,mp.prBndX,mp.prBndY,max(bctAll(:,9))+1);

nePrFxBnd=cell(size(fxBnd,1),1);
bndRng=[sort(fxBnd(:,[1,3]),2),sort(fxBnd(:,[2,4]),2)];

%% Check box crossing for each fixed boundary
for bnc=1:size(fxBnd,1)
    valCs=((sdRngCp(:,1)<bndRng(bnc,2)) & ...
        (sdRngCp(:,2)>bndRng(bnc,1)) & ...
        (sdRngCp(:,3)<bndRng(bnc,4)) & ...
        (sdRngCp(:,4)>bndRng(bnc,3)));
    neId=sdIdCp(valCs);
    nePrFxBnd{bnc}=[ones(numel(neId),1)*bnc,neId];
end

%% Create a pair of potential neighboring fixed boundary and cells
nePrFxBnd=nePrFxBnd(~cellfun('isempty',nePrFxBnd));
nePrFxBnd=cell2mat(nePrFxBnd);
nePrFxBnd=unique(nePrFxBnd,'rows');

end
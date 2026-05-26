%% Find all neighbor pairs, using disk radius

function nePr=bct_neiPairRectangleAll(bctAll,mp)

[~,~,~,sdRngCp,sdIdCp]=bct_seedCopy...
    (bctAll(:,1:2),bctAll(:,3),bctAll(:,9),bctAll(:,10:13),...
    mp.lx,mp.ly,mp.prBndX,mp.prBndY,max(bctAll(:,9))+1);

nePr=cell(mp.nFa,1);

for fac=1:mp.nFa
    valCs=((sdRngCp(:,1)<sdRngCp(fac,2)) & ...
        (sdRngCp(:,2)>sdRngCp(fac,1)) & ...
        (sdRngCp(:,3)<sdRngCp(fac,4)) & ...
        (sdRngCp(:,4)>sdRngCp(fac,3)));
    neId=sdIdCp(valCs);
    neId=neId(neId~=fac);
    nePr{fac}=sort([ones(numel(neId),1)*fac,neId],2);
end

nePr=nePr(~cellfun('isempty',nePr));
nePr=cell2mat(nePr);
nePr=unique(nePr,'rows');

end
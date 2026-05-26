%% Compute boundary edge length

function eLn=bct_boundaryEdgeLength(mvBnd,mp)

mvBndc=[mvBnd;mvBnd(1,:)];
tnVec=mvBndc(2:end,:)-mvBndc(1:end-1,:);

if mp.prBndX==1
    tnVec(:,1)=tnVec(:,1)-mp.lx.*floor(tnVec(:,1)./mp.lx+1/2);
end
if mp.prBndY==1
    tnVec(:,2)=tnVec(:,2)-mp.ly.*floor(tnVec(:,2)./mp.ly+1/2);
end

eLn=sqrt(sum(tnVec.^2,2));

end
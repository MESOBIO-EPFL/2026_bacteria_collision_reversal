%% Redistribute boundary vertices

function nmd=bct_movingBoundaryRedistribute(mvBnd,eLn,mp)

%% compute number of intermediate vertices
rfn=floor(sum(eLn)/mp.bctLn);
eLnNm=eLn/sum(eLn);

%% Get coordinates for vertices
mvBndc=[mvBnd;mvBnd(1,:)];
for vrc=2:size(mvBndc,1)
    if mp.prBndX==1
        mvBndc(vrc,1)=mvBndc(vrc,1)...
            -mp.lx.*floor((mvBndc(vrc,1)-mvBndc(vrc-1,1))./mp.lx+1/2);
    end
    if mp.prBndY==1
        mvBndc(vrc,2)=mvBndc(vrc,2)...
            -mp.ly.*floor((mvBndc(vrc,2)-mvBndc(vrc-1,2))./mp.ly+1/2);
    end
end

%% specify line parameter for piecewise linear function
erg=zeros(size(eLn,1)+1,1);
for ii=2:numel(erg)
    erg(ii)=sum(eLnNm(1:ii-1));
end

%% redistribute intermediate vertices
cnp=(1/(rfn+1))*((1:rfn+1)-1/2);
nmd=zeros(numel(cnp),2);
for ii=1:size(mvBndc,1)-1
    for jj=1:size(nmd,1)
        nmd(jj,:)=nmd(jj,:)+...
            (mvBndc(ii,:)+(cnp(jj)-erg(ii))/(erg(ii+1)-erg(ii))*...
            (mvBndc(ii+1,:)-mvBndc(ii,:)))*(heaviside(erg(ii+1)-cnp(jj))-...
            heaviside(erg(ii)-cnp(jj)));
    end
end    

if mp.prBndX==1
    nmd(:,1)=mod(nmd(:,1),mp.lx);
end
if mp.prBndY==1
    nmd(:,2)=mod(nmd(:,2),mp.ly);
end


end
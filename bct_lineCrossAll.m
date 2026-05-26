%% Compute whether two line segments intersect or not

function [cross,parVal]=bct_lineCrossAll(bctAll,nePr,mp)

nePrVec=(nePr.');
nePrVec=nePrVec(:);
nmx=numel(nePrVec);

%% locally adjust particle positions under periodic boundary condition
cnCrd=bctAll(nePrVec,1:2);
if mp.prBndX==1
    cnCrd(2:2:nmx,1)=cnCrd(2:2:nmx,1)-mp.lx.*...
        floor((cnCrd(2:2:nmx,1)-cnCrd(1:2:nmx,1))./mp.lx+1/2);
end
if mp.prBndY==1
    cnCrd(2:2:nmx,2)=cnCrd(2:2:nmx,2)-mp.ly.*...
        floor((cnCrd(2:2:nmx,2)-cnCrd(1:2:nmx,2))./mp.ly+1/2);
end

%% Compute x and y coordinates of two end points
xCrdIn=cnCrd(:,1)-bctAll(nePrVec,9)/2.*cos(bctAll(nePrVec,3));
xCrdFn=cnCrd(:,1)+bctAll(nePrVec,9)/2.*cos(bctAll(nePrVec,3));

yCrdIn=cnCrd(:,2)-bctAll(nePrVec,9)/2.*sin(bctAll(nePrVec,3));
yCrdFn=cnCrd(:,2)+bctAll(nePrVec,9)/2.*sin(bctAll(nePrVec,3));

%% Compute matrix coefficient for parameter calculation
cmXX=(xCrdFn(1:2:nmx)-xCrdIn(1:2:nmx));
cmXY=(xCrdIn(2:2:nmx)-xCrdFn(2:2:nmx));
cmYX=(yCrdFn(1:2:nmx)-yCrdIn(1:2:nmx));
cmYY=(yCrdIn(2:2:nmx)-yCrdFn(2:2:nmx));

detVal=cmXX.*cmYY-cmXY.*cmYX;

conX=(xCrdIn(2:2:nmx)-xCrdIn(1:2:nmx));
conY=(yCrdIn(2:2:nmx)-yCrdIn(1:2:nmx));

%% Compute parameter values
parVal(:,1)=(cmYY.*conX-cmXY.*conY)./detVal;
parVal(:,2)=(-cmYX.*conX+cmXX.*conY)./detVal;

%% Determine crossing or not
cross=(parVal(:,1)>0 & parVal(:,2)>0 & parVal(:,1)<1 & parVal(:,2)<1);

end
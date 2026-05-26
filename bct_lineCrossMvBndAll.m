%% Compute whether two line segments intersect or not

function [cross,parVal]=bct_lineCrossMvBndAll(bctAll,bndAll,nePrBnd,mp)

%% Get cell ID for crossing check
cnCrd=bctAll(nePrBnd(:,2),1:2);
bndAllc=[bndAll;bndAll(1,:)];

%% Compute x and y coordinates of two end points
xCrdIn=cnCrd(:,1)-bctAll(nePrBnd(:,2),9)/2.*cos(bctAll(nePrBnd(:,2),3));
xCrdFn=cnCrd(:,1)+bctAll(nePrBnd(:,2),9)/2.*cos(bctAll(nePrBnd(:,2),3));

xCrdBndIn=bndAllc(nePrBnd(:,1),1);
xCrdBndFn=bndAllc(nePrBnd(:,1)+1,1);

yCrdIn=cnCrd(:,2)-bctAll(nePrBnd(:,2),9)/2.*sin(bctAll(nePrBnd(:,2),3));
yCrdFn=cnCrd(:,2)+bctAll(nePrBnd(:,2),9)/2.*sin(bctAll(nePrBnd(:,2),3));

yCrdBndIn=bndAllc(nePrBnd(:,1),2);
yCrdBndFn=bndAllc(nePrBnd(:,1)+1,2);

% Correct local coordinate information
if mp.prBndX==1
    xCrdIn=xCrdIn-mp.lx.*floor((xCrdIn-xCrdBndIn)./mp.lx+1/2);
    xCrdFn=xCrdFn-mp.lx.*floor((xCrdFn-xCrdBndIn)./mp.lx+1/2);
    xCrdBndFn=xCrdBndFn-mp.lx.*floor((xCrdBndFn-xCrdBndIn)./mp.lx+1/2);
end
if mp.prBndY==1
    yCrdIn=yCrdIn-mp.ly.*floor((yCrdIn-yCrdBndIn)./mp.ly+1/2);
    yCrdFn=yCrdFn-mp.ly.*floor((yCrdFn-yCrdBndIn)./mp.ly+1/2);
    yCrdBndFn=yCrdBndFn-mp.ly.*floor((yCrdBndFn-yCrdBndIn)./mp.ly+1/2);
end

%% Compute matrix coefficient for parameter calculation
cmXX=(xCrdBndFn-xCrdBndIn);
cmXY=(xCrdIn-xCrdFn);
cmYX=(yCrdBndFn-yCrdBndIn);
cmYY=(yCrdIn-yCrdFn);

detVal=cmXX.*cmYY-cmXY.*cmYX;

conX=(xCrdIn-xCrdBndIn);
conY=(yCrdIn-yCrdBndIn);

%% Compute parameter values
parVal(:,1)=(cmYY.*conX-cmXY.*conY)./detVal;
parVal(:,2)=(-cmYX.*conX+cmXX.*conY)./detVal;

%% Determine crossing or not
cross=(parVal(:,1)>0 & parVal(:,2)>0 & parVal(:,1)<1 & parVal(:,2)<1);

end
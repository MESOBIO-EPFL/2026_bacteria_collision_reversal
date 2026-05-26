%% Compute whether two line segments intersect or not

function [cross,parVal]=bct_lineCrossFxBndAll(bctAll,bndAll,nePrBnd,mp)

%% Get cell ID for crossing check
cnCrd=bctAll(nePrBnd(:,2),1:2);

%% Compute x and y coordinates of two end points
xCrdIn=cnCrd(:,1)-bctAll(nePrBnd(:,2),9)/2.*cos(bctAll(nePrBnd(:,2),3));
xCrdFn=cnCrd(:,1)+bctAll(nePrBnd(:,2),9)/2.*cos(bctAll(nePrBnd(:,2),3));

xCrdBndIn=bndAll(nePrBnd(:,1),1);
xCrdBndFn=bndAll(nePrBnd(:,1),3);

yCrdIn=cnCrd(:,2)-bctAll(nePrBnd(:,2),9)/2.*sin(bctAll(nePrBnd(:,2),3));
yCrdFn=cnCrd(:,2)+bctAll(nePrBnd(:,2),9)/2.*sin(bctAll(nePrBnd(:,2),3));

yCrdBndIn=bndAll(nePrBnd(:,1),2);
yCrdBndFn=bndAll(nePrBnd(:,1),4);

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
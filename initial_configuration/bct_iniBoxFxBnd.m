%% Define boundary edges

function fxBnd=bct_iniBoxFxBnd(mp,exWth)

bndVr=[(mp.lx-exWth)/2,mp.ly+10;...
    (mp.lx-exWth)/2,mp.ly;0,mp.ly;...
    0,0;mp.lx,0;mp.lx,mp.ly;...
    (mp.lx+exWth)/2,mp.ly;...
    (mp.lx+exWth)/2,mp.ly+10];

fxBnd=zeros(size(bndVr,1)-1,6);
fxBnd(:,1:2)=bndVr(1:end-1,:);
fxBnd(:,3:4)=bndVr(2:end,:);

tnVec=fxBnd(:,3:4)-fxBnd(:,1:2);
nmVec=[-tnVec(:,2),tnVec(:,1)];
nmVecMag=sqrt(sum(nmVec.^2,2));
nmVec=nmVec./nmVecMag;
fxBnd(:,5:6)=nmVec;

end
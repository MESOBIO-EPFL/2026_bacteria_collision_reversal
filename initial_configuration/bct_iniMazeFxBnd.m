%% Define initial boundary condition

function fxBnd=bct_iniMazeFxBnd(mp)

%% Define fixed boundary on the left end
bndVr=[47.5-mp.mzLn*sin(mp.mzAng),72.5+mp.mzLn*cos(mp.mzAng);...
    47.5,72.5;25,72.5;25,60;82.5,60;82.5,42.5;12.5,42.5;...
    12.5,45;80,45;80,57.5;22.5,57.5;22.5,72.5;0,72.5;0,0;...
    32.5,0;32.5,12.5;12.5,12.5;12.5,15;35,15;35,0;...
    95,0;95,12.5;80,12.5;80,27.5;17.5,27.5;17.5,30;...
    82.5,30;82.5,15;95,15;95,72.5;60,72.5;...
    60+mp.mzLn*sin(mp.mzAng),72.5+mp.mzLn*cos(mp.mzAng)];

fxBnd=zeros(size(bndVr,1)-1,6);
fxBnd(:,1:2)=bndVr(1:end-1,:);
fxBnd(:,3:4)=bndVr(2:end,:);

tnVec=fxBnd(:,3:4)-fxBnd(:,1:2);
nmVec=[-tnVec(:,2),tnVec(:,1)];
nmVecMag=sqrt(sum(nmVec.^2,2));
nmVec=nmVec./nmVecMag;
fxBnd(:,5:6)=nmVec;

end
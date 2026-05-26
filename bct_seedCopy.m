%% Copy original seed to surrounding boxes
function [sdCp,sdOrnCp,sdLenCp,sdRngCp,sdIdCp]=...
    bct_seedCopy(sd,sdOrn,sdLen,sdRng,lx,ly,prBndX,prBndY,rd)

%% Define copy condition for different boundary condition
if prBndX==1 && prBndY==1
    cpBx=[0,0;lx,ly;lx,0;lx,-ly;0,ly;0,-ly;-lx,ly;-lx,0;-lx,-ly];
elseif prBndX==1 && prBndY~=1
    cpBx=[0,0;lx,0;-lx,0];
elseif prBndX~=1 && prBndY==1
    cpBx=[0,0;0,ly;0,-ly];
else
    cpBx=[0,0];
end

%% Copy seed to adjacent boxes
sdCnt=size(sd,1);
sdCp=zeros(size(cpBx,1)*sdCnt,2);
sdOrnCp=repmat(sdOrn,size(cpBx,1),1);
sdLenCp=repmat(sdLen,size(cpBx,1),1);
sdIdCp=repmat((1:sdCnt).',size(cpBx,1),1);
sdRngCp=zeros(size(cpBx,1)*sdCnt,4);

for cbc=1:size(cpBx,1)
    sdCp(1+(cbc-1)*sdCnt:cbc*sdCnt,:)=sd+repmat(cpBx(cbc,:),sdCnt,1);  
    sdRngCp(1+(cbc-1)*sdCnt:cbc*sdCnt,:)=...
        [sdRng(:,1:2)+cpBx(cbc,1),sdRng(:,3:4)+cpBx(cbc,2)];
end

%% Take cells only within a margin of bndTh
bndTh=1.1*rd;

if prBndX==1 && prBndY==1
    valCs=(sdCp(:,1)>-bndTh & sdCp(:,1)<lx+bndTh & ...
        sdCp(:,2)>-bndTh & sdCp(:,2)<ly+bndTh);
elseif prBndX==1 && prBndY==0
    valCs=(sdCp(:,1)>-bndTh & sdCp(:,1)<lx+bndTh);
elseif prBndX==0 && prBndY==1
    valCs=(sdCp(:,2)>-bndTh & sdCp(:,2)<ly+bndTh);
else
    valCs=true(sdCnt,1);
end

sdCp=sdCp(valCs,:);
sdOrnCp=sdOrnCp(valCs);
sdLenCp=sdLenCp(valCs);
sdIdCp=sdIdCp(valCs);
sdRngCp=sdRngCp(valCs,:);

end
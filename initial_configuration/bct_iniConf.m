%% Create initial configuration (Periodic)

function [bctAll,dtAll,mpN]=bct_iniConf(mp)

%% Create output variables
bctAll=zeros(mp.nFa,13);
dtAll=cell(floor(mp.ttTm/mp.tmSt/mp.svPt),1);
mpN=mp;

%% Randomly assign cell length from bctLn to 2*bctLn+1
bctAll(:,9)=mpN.bctLn+(mpN.bctLn+1)*rand(mpN.nFa,1);
mpN.lx=sqrt(mpN.nFa/mpN.vf/mpN.rto*mean(bctAll(:,9)+pi/4));
mpN.ly=mpN.rto*mpN.lx;

if mpN.tsRev<=0
    bctAll(:,6)=0;
else
    bctAll(:,6)=rand(mp.nFa,1)*mpN.tsRev;
end

%% Randomly generate initial position and orientation
bctAll(:,1)=rand(mpN.nFa,1)*mpN.lx;
bctAll(:,2)=rand(mpN.nFa,1)*mpN.ly;
bctAll(:,3)=rand(mpN.nFa,1)*2*pi;
bctAll(:,4)=ones(mpN.nFa,1);

%% Compute cell range 
bctAll(:,10:13)=bct_cellRangeCompute...
    (bctAll(:,1:2),bctAll(:,3),bctAll(:,9),mpN.rngBf);
mpN.sdRngRef=bctAll(:,10:13);

end
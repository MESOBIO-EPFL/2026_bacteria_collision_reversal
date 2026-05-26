%% Create initial configuration (Periodic)

function [bctAll,dtAll,mpN]=bct_iniConfMaze(mp)

%% Creat initial configuration
bctAll=zeros(mp.nFa,13);
if mp.tsRev<=0
    bctAll(:,6)=zeros(mp.nFa,1)*mp.tsRev;
else
    bctAll(:,6)=rand(mp.nFa,1)*mp.tsRev;
end

% Randomly generate initial position and orientation
bctAll(:,1)=(1.5*mp.bctLn+35)+(45-3*mp.bctLn)*rand(mp.nFa,1);
bctAll(:,2)=1.5*mp.bctLn+rand(mp.nFa,1)*(27.5-3*mp.bctLn);

bctAll(:,3)=rand(mp.nFa,1)*2*pi;
bctAll(:,4)=ones(mp.nFa,1);

dtAll=cell(floor(mp.ttTm/mp.tmSt/mp.svPt),1);

% Randomly assign cell length from bctLn to 2*bctLn+1
bctAll(:,9)=mp.bctLn+(mp.bctLn+1)*rand(mp.nFa,1);

% Compute cell range 
bctAll(:,10:13)=bct_cellRangeCompute...
    (bctAll(:,1:2),bctAll(:,3),bctAll(:,9),mp.rngBf);

% Save reference cell range in mp file
mpN=mp;
mpN.sdRngRef=bctAll(:,10:13);

end
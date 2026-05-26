%% Cell division protocol

function [bctAllN,mpN,dvOcc]=bct_cellDivision(bctAll,mp)

% Check whether any cells are longer than the division criterion
mpN=mp;
bctAllN=bctAll;

dvChk=find(bctAll(:,9)>2*mp.bctLn+1);
dvOcc=0;

if ~isempty(dvChk)
    dvOcc=1;

    dvmx=numel(dvChk);
    newCell=zeros(dvmx,13);

    %% Update cell position
    bctAllN(dvChk,1)=bctAllN(dvChk,1)-(bctAllN(dvChk,9)/4+1/4)...
        .*cos(bctAllN(dvChk,3));
    bctAllN(dvChk,2)=bctAllN(dvChk,2)-(bctAllN(dvChk,9)/4+1/4)...
        .*sin(bctAllN(dvChk,3));

    newCell(:,1)=bctAll(dvChk,1)+(bctAll(dvChk,9)/4+1/4)...
        .*cos(bctAll(dvChk,3));
    newCell(:,2)=bctAll(dvChk,2)+(bctAll(dvChk,9)/4+1/4)...
        .*sin(bctAll(dvChk,3));

    %% Update orientation, polarity, moment, reversal time and count.
    newCell(:,3:8)=bctAllN(dvChk,3:8);

    if mp.tsRev>0
        newCell(:,6)=newCell(:,6)+0.02*rand(dvmx,1)*mp.tsRev;
    end

    %% Update cell length    
    bctAllN(dvChk,9)=bctAllN(dvChk,9)/2-1/2;
    newCell(:,9)=bctAllN(dvChk,9);

    %% Merge new cell array into existing cell matrix
    bctAllN=[bctAllN;newCell];

    %% Compute cell range 
    bctAllN(:,10:13)=bct_cellRangeCompute(bctAllN(:,1:2),...
        bctAllN(:,3),bctAllN(:,9),mp.rngBf);
    
    %% Update model parameters
    mpN.nFa=mpN.nFa+dvmx;
    mpN.sdRngRef=bctAllN(:,10:13);
    
    mpN.mtCellId=[mpN.mtCellId;mpN.mtCellId(dvChk)];
    mpN.tmSclFrt=[mpN.tmSclFrt;mpN.tmSclFrt(dvChk)];
    mpN.dvTm=[mpN.dvTm;(0.8+0.4*rand(numel(dvChk),1))*mp.tsDiv];

    if sum(bctAllN(:,9)+pi/4)/mp.lx/mp.ly>=mp.vfMx
        mpN.dvOn=0;
    end
end

end
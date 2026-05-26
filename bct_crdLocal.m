%% Correct vertex coordinates. 
% If vertex crosses the boundary, correct it. 

function vrF=bct_crdLocal(vrL,lx,ly,prBndX,prBndY)

%% final vertex coordinats
vrF=vrL;

if prBndX==1
    vrF(:,1)=vrF(:,1)-lx.*floor((vrF(:,1)-vrF(1,1))./lx+1/2);
end

if prBndY==1
    vrF(:,2)=vrF(:,2)-ly.*floor((vrF(:,2)-vrF(1,2))./ly+1/2);
end

end
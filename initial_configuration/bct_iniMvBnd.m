%% Define initial boundary condition

function mvBnd=bct_iniMvBnd(mp)

%% Define moving boundary on the right end
vrCnt=floor(mp.ly/mp.bctLn);
mvBnd=[ones(vrCnt,1)*(mp.lx+2*mp.bctLn+3),...
    (0:mp.ly/vrCnt:mp.ly*(vrCnt-1)/vrCnt+0.01).'];
end
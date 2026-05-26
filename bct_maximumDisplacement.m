%% Check maximum displacement in both x and y for all cells

function mxDis=bct_maximumDisplacement(bctAll,mp)

dis=mp.sdRngRef-bctAll(:,10:13);

if mp.prBndX==1
    dis(:,1:2)=dis(:,1:2)-mp.lx.*floor((dis(:,1:2))./mp.lx+1/2);
end

if mp.prBndY==1
    dis(:,3:4)=dis(:,3:4)-mp.ly.*floor((dis(:,3:4))./mp.ly+1/2);
end

dis=abs(dis);
mxDis=max(dis(:));

end
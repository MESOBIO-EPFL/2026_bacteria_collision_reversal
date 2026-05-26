%% Compute cell range for a given position and orientation

function cellRange=bct_cellRangeCompute(pos,orien,len,rngBf)

cellRange=[pos(:,1)-len/2.*abs(cos(orien))-1/2-rngBf,...
    pos(:,1)+len/2.*abs(cos(orien))+1/2+rngBf,...
    pos(:,2)-len/2.*abs(sin(orien))-1/2-rngBf,...
    pos(:,2)+len/2.*abs(sin(orien))+1/2+rngBf];

end
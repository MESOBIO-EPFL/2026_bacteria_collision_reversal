%% Save data

function bct_saveDataFxMvBnd(dtAll,dtBndAll,fxBnd,mp,csNm,odr,rpc)

%% Create output directory
if ~exist(odr,'dir')
    mkdir(odr);
end

%% Save data in Mat file format
save(sprintf('%s%s_r%d_data.mat',odr,csNm,rpc),'dtAll');
save(sprintf('%s%s_r%d_mvBnd_data.mat',odr,csNm,rpc),'dtBndAll');
save(sprintf('%s%s_r%d_fxBnd_data.mat',odr,csNm,rpc),'fxBnd');
save(sprintf('%s%s_r%d_mp.mat',odr,csNm,rpc),'mp');


end
%% Save data

function bct_saveDataFxBnd(dtAll,fxBnd,mp,csNm,odr,rpc)

%% Create output directory
% odr='./data/';
if ~exist(odr,'dir')
    mkdir(odr);
end

%% Save data in Mat file format
save(sprintf('%s%s_r%d_data.mat',odr,csNm,rpc),'dtAll');
save(sprintf('%s%s_r%d_bnd_data.mat',odr,csNm,rpc),'fxBnd');
save(sprintf('%s%s_r%d_mp.mat',odr,csNm,rpc),'mp');


end
%% Save data

function bct_saveData(dtAll,odr,mp,csNm,rpc)

%% Create output directory
if ~exist(odr,'dir')
    mkdir(odr);
end

%% Save data in Mat file format
save(sprintf('%s%s_r%d_data.mat',odr,csNm,rpc),'dtAll');
save(sprintf('%s%s_r%d_mp.mat',odr,csNm,rpc),'mp');


end
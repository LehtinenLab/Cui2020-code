fdir = 'C:\Users\LehtinenLab\Dropbox\Jin\20200403_multi_analysis';

CCL2name = {
    'EmbryonicCX3CR1_190821_1_4Dchunk',
    'EmbryonicCX3CR1_190821_2_4Dchunk',
    'EmbryonicCX3CR1_190814_3_4Dchunk',
    'EmbryonicCX3CR1_190904_2_4Dchunk',
    'EmbryonicCX3CR1_200715_2_4Dchunk'
};

CREname = {
    'EmbryonicCX3CR1_190925_1_4Dchunk',
    'EmbryonicCX3CR1_190703_4_4Dchunk',
    'EmbryonicCX3CR1_190911_3_4Dchunk',
    'EmbryonicCX3CR1_200715_1_4Dchunk'
};

%%
CCL2_sorted_all = [];
for i = 1:numel(CCL2name)
    CCL2data = load(strcat(fdir,filesep,CCL2name{i},'.data'),'data','-mat');
    CCL2data = CCL2data.data;
    if strcmp(class(CCL2data),'table')
        CCL2data = table2struct(CCL2data);
    end
    [CCL2_sorted_temp,CCL2_M{i}] = MakeSortedMotilityMap(CCL2data);
    cellnum(i) = size(CCL2_sorted_temp,1);
    CCL2_sorted_all = cat(1,CCL2_sorted_all,CCL2_sorted_temp);
end
cellnumcum = cumsum(cellnum);

clo = min(CCL2_sorted_all(:));
chi = max(CCL2_sorted_all(:));

CCL2_sorted_M = [];
for i = 1:numel(CCL2name)
    CCL2_sorted_M = cat(2,CCL2_sorted_M,rescale(CCL2_M{i},clo,chi));
end
% CCL2_sorted_M = rescale(CCL2_sorted_M,clo,chi);



figure,
subplot(1,2,2)
imagesc(CCL2_sorted_all,[clo,chi]);
% imagesc(cat(2,CCL2_sorted_all,CCL2_sorted_M'),[clo,chi]);
colormap(gray);
colorbar;
title('CCL2 motility');
hold on
for i = 1:numel(cellnum)-1
    plot([0,size(CCL2_sorted_all,2)+0.5],[cellnumcum(i)+0.5,cellnumcum(i)+0.5],'--r','LineWidth',2);
end
% plot([size(CCL2_sorted_all,2)+0.5,size(CCL2_sorted_all,2)+0.5],[0,size(CCL2_sorted_all,1)],'--g','LineWidth',2);
hold off

CRE_sorted_all = [];

for i = 1:numel(CREname)
    CREdata = load(strcat(fdir,filesep,CREname{i},'.data'),'data','-mat');
    CREdata = CREdata.data;
    if strcmp(class(CREdata),'table')
        CREdata = table2struct(CREdata);
    end
    [CRE_sorted_temp, CRE_M{i}] = MakeSortedMotilityMap(CREdata);
    cellnum(i) = size(CRE_sorted_temp,1);
    CRE_sorted_all = cat(1,CRE_sorted_all,CRE_sorted_temp);
end
cellnumcum = cumsum(cellnum);

CRE_sorted_M = [];
for i = 1:numel(CRE_M)
    CRE_sorted_M = cat(2,CRE_sorted_M,rescale(CRE_M{i},clo,chi));
end


subplot(1,2,1),
imagesc(CRE_sorted_all,[clo,chi]);
% imagesc(cat(2,CRE_sorted_all,CRE_sorted_M'),[clo,chi]);
colormap(gray);
colorbar;
title('CRE motility');
hold on
for i = 1:numel(cellnum)-1
    plot([0,size(CRE_sorted_all,2)+0.5],[cellnumcum(i)+0.5,cellnumcum(i)+0.5],'--r','LineWidth',2);
end
% plot([size(CRE_sorted_all,2)+0.5,size(CRE_sorted_all,2)+0.5],[0,size(CRE_sorted_all,1)],'--g','LineWidth',2);
hold off


function [CCL2_sorted,M] = MakeSortedMotilityMap(BPdata)
    A = [BPdata(:).cell_bin];
    A = reshape(A,size(BPdata(1).cell_bin,1),size(BPdata(1).cell_bin,2),[]);
    A_sz = squeeze(sum(sum(A,1),2));
    idx = A_sz > 50;
    
    ii = find(idx);
    for i = 1:numel(ii)
        M_temp = BPdata(ii(i)).cell_DFT;
        M(i) = mean(M_temp(:));
    end
    BP_std_norm = [BPdata(idx).BP_ring_norm];
    BP_std_norm = reshape(BP_std_norm,[],sum(idx));
    BP_std_norm = BP_std_norm';
    pks = max(BP_std_norm,[],2);
    [~,sort_pks] = sort(pks);
    sort_pks = flip(sort_pks);
    M = M(sort_pks);
    CCL2_sorted = BP_std_norm(sort_pks,:);
end
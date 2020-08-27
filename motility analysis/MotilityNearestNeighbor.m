fdir = 'C:\Users\LehtinenLab\Dropbox\Jin\20200403_multi_analysis';


fnamelist = {
    'EmbryonicCX3CR1_190821_1_4Dchunk',
    'EmbryonicCX3CR1_190821_2_4Dchunk',
    'EmbryonicCX3CR1_190814_3_4Dchunk',
    'EmbryonicCX3CR1_190904_2_4Dchunk',
    'EmbryonicCX3CR1_200715_2_4Dchunk',
    'EmbryonicCX3CR1_190925_1_4Dchunk',
    'EmbryonicCX3CR1_190703_4_4Dchunk',
    'EmbryonicCX3CR1_190911_3_4Dchunk',
    'EmbryonicCX3CR1_200715_1_4Dchunk'
};

C = {'CCL2','CCL2','CCL2','CCL2','CCL2','CRE','CRE','CRE','CRE'};

scale = [0.5,0.5,5];
%% calculate NN distances and save as table
for i = 1:numel(fnamelist)
% for i = 5
    fname = fnamelist{i};
    data = load(strcat(fdir,filesep,fname,'.data'),'data','-mat');
    try
        data = struct2table(data.data);
        centroids = data.centroid;
        ncells = size(data,1);
        NN_dist_um = zeros(1,ncells);
        for j = 1:ncells
            c = repmat(centroids(j,:),[ncells,1]);
            d = centroids - c;
            d = d.* repmat(scale,[ncells,1]);
            dist = vecnorm(d,2,2);
            NN_dist_um(j) = min(dist(dist>0));
        end
        data = addvars(data,NN_dist_um','NewVariableNames','NN_dist_um');
        save(strcat(fdir,filesep,fnamelist{i},'.data'),'data');
        fprintf('saved file %d of %d\n',i,numel(fnamelist));
    catch
        fprintf('file %d of %d already converted\n',i,numel(fnamelist));
    end
end

%%
w = waitbar(0,'building data table...');
for i = 1:numel(fnamelist)
% for i = 1
    fname = fnamelist{i};
    data = load(strcat(fdir,filesep,fname,'.data'),'data','-mat');
    data = data.data;
    
    cell_idx = data.idx;

%     displacement = data.disp_total*scale(1);
    R_disp = [data.R_disp{:}];
    C_disp = [data.C_disp{:}];
    displacement_step = scale(1)*sqrt(diff(R_disp,1,1).^2 + diff(C_disp,1,1).^2);
    displacement = sum(displacement_step,1)';

    NN = data.NN_dist_um;
    brightness = data.cell_brightness;
    motility = data.BP_ring_norm;
    condition = repmat(C(i),size(cell_idx));
    [~,condition_idx] = ismember(C(i),unique(C));
    condition_idx = repmat(condition_idx,size(cell_idx));
    sample_idx = repmat(i,size(cell_idx));
    
    
    T_temp = table(condition,condition_idx,sample_idx,cell_idx,displacement,NN,brightness,motility);
    try
        T = vertcat(T,T_temp);
    catch
        T = T_temp;
    end
    waitbar(i/numel(fnamelist));
end
close(w);
%%
motility_mean = nanmean(T.motility,1);
[~,motility_pk_loc] = max(motility_mean);
T.motility_pk = T.motility(:,motility_pk_loc);
%%
c1 = T.condition;
c1_L = unique(c1);
c2 = T.sample_idx;
labels = unique(strcat(c1," ",num2str(c2)));
idx_1 = ismember(c1,c1_L(1));
idx_2 = ismember(c1,c1_L(2));

figure,

subplot(2,3,1);
% gscatter(T.NN,T.displacement,T.condition,'rbgk');

gscatter(T.NN(idx_1),T.displacement(idx_1),c2(idx_1),'r','+*xsd^v',[],'off');
hold on
gscatter(T.NN(idx_2),T.displacement(idx_2),c2(idx_2),'b','+*xsd^v',[],'off');
legend(labels);
xlabel('Nearest-neighbor distance');
ylabel('Cell displacement');


subplot(2,3,2);
% gscatter(T.brightness,T.displacement,T.condition,'rbgk');
gscatter(T.brightness(idx_1),T.displacement(idx_1),c2(idx_1),'r','+*xsd^v',[],'off');
hold on
gscatter(T.brightness(idx_2),T.displacement(idx_2),c2(idx_2),'b','+*xsd^v',[],'off');
legend(labels);

xlabel('Cell body brightness');
ylabel('Cell displacement');


subplot(2,3,3);
% gscatter(T.NN,T.brightness,T.condition,'rbgk');
gscatter(T.NN(idx_1),T.brightness(idx_1),c2(idx_1),'r','+*xsd^v',[],'off');
hold on
gscatter(T.NN(idx_2),T.brightness(idx_2),c2(idx_2),'b','+*xsd^v',[],'off');
legend(labels);
xlabel('Nearest-neighbor distance');
ylabel('Cell body brightness');

subplot(2,3,4);
% gscatter(T.NN,T.motility_pk,T.condition,'rbgk');
gscatter(T.NN(idx_1),T.motility_pk(idx_1),c2(idx_1),'r','+*xsd^v',[],'off');
hold on
gscatter(T.NN(idx_2),T.motility_pk(idx_2),c2(idx_2),'b','+*xsd^v',[],'off');
legend(labels);
xlabel('Nearest-neighbor distance');
ylabel('Peak motility');

subplot(2,3,5);
% gscatter(T.brightness,T.motility_pk,T.condition,'rbgk');
gscatter(T.brightness(idx_1),T.motility_pk(idx_1),c2(idx_1),'r','+*xsd^v',[],'off');
hold on
gscatter(T.brightness(idx_2),T.motility_pk(idx_2),c2(idx_2),'b','+*xsd^v',[],'off');
legend(labels);
xlabel('Cell body brightness');
ylabel('Peak motility');

subplot(2,3,6);
% gscatter(T.displacement,T.motility_pk,T.condition,'rbgk');
gscatter(T.displacement(idx_1),T.motility_pk(idx_1),c2(idx_1),'r','+*xsd^v',[],'off');
hold on
gscatter(T.displacement(idx_2),T.motility_pk(idx_2),c2(idx_2),'b','+*xsd^v',[],'off');
legend(labels);
xlabel('Cell Displacement');
ylabel('Peak motility');

%% log-log plot of displacement vs motility
figure,
gscatter(T.displacement(idx_1),T.motility_pk(idx_1),c2(idx_1),'r','+*xsd^v',[],'off');
hold on
gscatter(T.displacement(idx_2),T.motility_pk(idx_2),c2(idx_2),'b','+*xsd^v',[],'off');
set(gca,'xscale','log')
set(gca,'yscale','log')
ylim([0.045,0.4]);
legend(labels,'Location','northeastoutside');
xlabel('Cell displacement (um)');
ylabel('Peak motility');

%% displacement/motility grouped by FOV

for i = 1:max(T.sample_idx)
    mean_disp(i) = mean(T.displacement(T.sample_idx == i));
    mean_motility(i) = mean(T.motility(T.sample_idx == i));
    mean_condition(i) = mean(T.condition_idx(T.sample_idx == i));
end

figure,

subplot(1,2,1);
gscatter(mean_disp,mean_motility,mean_condition,'rb');
xlabel('Cell Displacement');
ylabel('Peak motility');
axis square
legend('CCL2','CRE','Location','northoutside');

subplot(1,2,2);
gscatter(mean_disp,mean_motility,mean_condition,'rb');
set(gca,'xscale','log')
set(gca,'yscale','log')
xlabel('Cell Displacement');
ylabel('Peak motility');
axis square
legend('CCL2','CRE','Location','northoutside');

sgtitle('displacement vs. motility (mean FOV)');
javaaddpath 'C:\Program Files\MATLAB\R2018a\java\mij.jar'
javaaddpath 'C:\Program Files\MATLAB\R2018a\java\ij-1.52a.jar'
javaaddpath 'C:\Users\LehtinenLab\Dropbox\AndermannLab\users\Fred\TurboRegHL_.jar'

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

for i = 1:numel(fnamelist)
    fname = fnamelist{i};
    data = load(strcat(fdir,filesep,fname,'.data'),'data','-mat');
    data = data.data;
    %%
    
    Nc = size(data,1);
    Nx = size(data.cell_DFT{1},1);
    Ny = size(data.cell_DFT{1},2);
    Nt = size(data.cell_DFT{1},3);

%     frame = zeros(ceil(Nc/4)*(Ny+5),4*(3*Nx+10),3,Nt);
%     xpos = repmat(1:3*Nx+10:size(frame,2),[1,ceil(Nc/4)]);
%     ypos = repmat(1:Ny+5:size(frame,1),[1,4]);
%     xpos = sort(xpos);

    frame = zeros(10*(Ny+5),ceil(Nc/10)*(3*Nx+10),3,Nt);
    ypos = repmat(1:Ny+5:size(frame,1),[1,ceil(Nc/10)]);
    xpos = repmat(1:3*Nx+10:size(frame,2),[1,10]);
    xpos = sort(xpos);


    
    motility = mean(data.BP_ring_norm,2);
    [~,motility_order] = sort(motility,'descend');
    %
    for j = 1:Nc
%     for j = 12
        DFT_im = data.cell_DFT_unmask{motility_order(j)};
        DFT_rgb = repmat(rescale(DFT_im),[1,1,1,3]);
        DFT_rgb = permute(DFT_rgb,[1,2,4,3]);
        
        BP_im = data.BP_filt_unmask{motility_order(j)};
        BP_rgb = repmat(rescale(BP_im),[1,1,1,3]);
        BP_rgb = permute(BP_rgb,[1,2,4,3]);

        % make a movie of the index number
        num = zeros(size(BP_rgb,1),size(BP_rgb,2));
        F = figure;
        imshow(num);
        hold on;
        t=text(Nx/2,Ny/2,num2str(j),'HorizontalAlignment','center');
        t.Color = 'white';
        t.FontSize = 20;
        numim = getframe(F);
        set(gca,'position',[0 0 1 1],'units','normalized')
        close(F);
        
        numim = numim.cdata;
        numim = imresize3(numim,[Nx,Ny,3]);
        numim = repmat(numim,[1,1,1,Nt]);
        
        %concatenate the three movies
        im = cat(2,rescale(numim),rescale(DFT_rgb),rescale(BP_rgb));
        frame(ypos(j):ypos(j)+size(im,1)-1,xpos(j):xpos(j)+size(im,2)-1,:,:) = im;
    end
%     implay(frame);
    
    montage_name = strcat(fdir,filesep,fname,'_montage_20200720.tif');
    pipe.io.writeTiff(frame,montage_name);
end
%%




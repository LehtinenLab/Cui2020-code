%% Initilize
javaaddpath 'C:\Program Files\MATLAB\R2018a\java\mij.jar'
javaaddpath 'C:\Program Files\MATLAB\R2018a\java\ij-1.52a.jar'
javaaddpath 'C:\Users\LehtinenLab\Dropbox\AndermannLab\users\Fred\TurboRegHL_.jar'

%put in some identifying informationm
mouse = 'EmbryonicCX3CR1';
date = '200715'; %YYMMDD format
run = 2;
ftype = 'sbxz';
server = 'Zombie';

outdir = 'C:\Users\LehtinenLab\Dropbox\Jin\20200403_multi_analysis';

%size in microns
volsize = [255,255,400];

%desired scale (microns/px)
scale = [0.497,0.497,5];

%total time
T_total = 2537.389;
%%
fdir = pipe.lab.datedir(mouse,date,server);
path = pipe.lab.datapath(mouse,date,run,ftype,server);

info = pipe.io.sbxInfo(path);
Nchan = info.nchan;
Nx = info.sz(1);
Ny = info.sz(2);
Nz = info.otlevels;
Nvol = floor(info.nframes / info.otlevels);
%NUMBER OF FRAMES FOR 5 MINUTES
Nt = round(5*60*(Nvol/T_total));

[Zproj,~,~,~] = zproj_reg(path,1,25,-1,1:Nz,'mtype','.sbxz','proj_type','mean');

meanproj = squeeze(mean(Zproj,4));
Fmask = figure;
imshowpair(squeeze(meanproj(1,:,:)),squeeze(meanproj(2,:,:)),'ColorChannels',[1,2,0]);
roi = drawpolygon;
mask = uint16(createMask(roi));
b = regionprops(mask,'BoundingBox');
b = round(b.BoundingBox);
close(Fmask);

[~,RS,CS,~] = zproj_reg(path,1,Nt,-1,1:Nz,'mtype','.sbxz','proj_type','mean');

tic


mov = pipe.io.sbxRead(path,1,Nz*Nt,-1,[]);
mov = reshape(mov,Nchan,Nx,Ny,[],Nt);
mov = mov(:,b(2):b(2)+b(4),b(1):b(1)+b(3),:,:);
for j = 1:Nz
    for c = 1:Nchan
        for t = 1:Nt
            im = squeeze(mov(c,:,:,j,t));
            im = im.*mask(b(2):b(2)+b(4),b(1):b(1)+b(3));
            mov(c,:,:,j,t) = imtranslate(im,[CS(t),RS(t)]);
        end
    end
end
%resize such that each pixel is 0.5um, and each zplane is 5microns
mov = imresizen(double(mov),[1,volsize(1)/Nx/scale(1),volsize(2)/Ny/scale(2),volsize(3)/(Nz-1)/scale(3),1]);
mov = uint16(mov(2,:,:,:,:));

%     [~,~,M2] = EmbryonicMaxMedRatio(mov);
write2chanTiff(mov,strcat(outdir,filesep,mouse,'_',num2str(date),'_',num2str(run),'_','4Dchunk.tif'));

clear mov;

toc;

% save(strcat(fdir,filesep,'MaxMed_ratio.mat'),'R','G');
javaaddpath 'C:\Program Files\MATLAB\R2018a\java\mij.jar'
javaaddpath 'C:\Program Files\MATLAB\R2018a\java\ij-1.52a.jar'
javaaddpath 'C:\Users\LehtinenLab\Dropbox\AndermannLab\users\Fred\TurboRegHL_.jar'

%put in some identifying informationm
mouse = '918-122-324';
date = '200507'; %YYMMDD format
run = 2;
ftype = 'sbxz';
server = 'Zombie';
path = pipe.lab.datapath(mouse,date,run,ftype,server);

info = pipe.io.sbxInfo(path);
Nchan = info.nchan;
Nx = info.sz(1);
Ny = info.sz(2);
Nz = info.otlevels;
Nt = floor(info.nframes / info.otlevels);

t1 = 1; %T start
t2 = 2000; %T end

%show average of 10 volumes after t1
startvol = pipe.imread(path,t1*Nz+1,10*Nz,-1,[]);

refchan = 2;
%% draw ROI rectangle
meanproj = mean(startvol,4);
F = figure;
imshowpair(squeeze(meanproj(1,:,:)),squeeze(meanproj(2,:,:)),'ColorChannels',[1,2,0]);
roi = drawrectangle;
Rect = round(roi.Position);
close(F);
%% show overlaid region to check
ROIoverlay_rgb = zeros(Nx,Ny,3);
ROIoverlay_rgb(:,:,1) = rescale(squeeze(meanproj(1,:,:)));
ROIoverlay_rgb(:,:,2) = rescale(squeeze(meanproj(2,:,:)));
ROIoverlay_rgb(Rect(2):Rect(2)+Rect(4),Rect(1):Rect(1)+Rect(3),3) = ones(Rect(4)+1,Rect(3)+1);
figure,imshow(ROIoverlay_rgb);
%% See average z stack to set Z1 and Z2 below
V_whole = reshape(startvol,Nchan,Nx,Ny,Nz,[]);
V = mean(V_whole(:,Rect(2):Rect(2)+Rect(4),Rect(1):Rect(1)+Rect(3),:,:),5);
implay2chan(V,[]);
%% Load small section

%define z start and end points based on video above
z1 = 3; %Z start
z2 = 42; %Z end

clear startvol V_whole;
cropmov = zeros(Nchan, Rect(4)+1, Rect(3)+1, z2-z1+1, t2-t1+1);
w = waitbar(0,'loading planes...');
for z = z1:z2
    tempslice = pipe.imread(path,t1,t2-t1+1,-1,z);
    cropmov(:,:,:,z-z1+1,:) = tempslice(:,Rect(2):Rect(2)+Rect(4),Rect(1):Rect(1)+Rect(3),:);
    waitbar((z-z1)/(z2-z1));
end
close(w);
%
cropmov_proj = squeeze(max(cropmov,[],4));

%% Fine tune registration for just this area
[~,R,C,~] = zproj_reg(path, t1, t2-t1+1, -1, z1:z2,'mtype',strcat('.',ftype),'proj_type','mean','refchan',refchan);

%% load and apply fine tuned registration

cropmov = zeros(Nchan, Rect(4)+1, Rect(3)+1, z2-z1+1, t2-t1+1);
pathz = pipe.lab.datapath(mouse,date,run,ftype,server);
w = waitbar(0,'loading planes...');
for z = z1:z2
    tempslice = pipe.imread(pathz,t1,t2-t1+1,-1,z);
    slice_reg = zeros(size(tempslice));
    for c = 1:2
        for t = 1:numel(R)
            slice_reg(c,:,:,t) = imtranslate(squeeze(tempslice(c,:,:,t)),[C(t),R(t)]);
        end
    end
    cropmov(:,:,:,z-z1+1,:) = slice_reg(:,Rect(2):Rect(2)+Rect(4),Rect(1):Rect(1)+Rect(3),:);
    waitbar((z-z1)/(z2-z1));
end
close(w);

%% save registered 4D region as tif
savepath = strcat(pipe.lab.rundir(mouse,date,run,server),'_4D_section.tif');
write2chanTiff(uint16(cropmov),savepath);

%% xy - zy montage -- Shows montage of long axis 

xy_proj = squeeze(mean(cropmov,4));
zy_proj = squeeze(mean(cropmov,3));
zy_proj = flip(zy_proj,3);

newsize = [size(zy_proj,1),size(zy_proj,2),2*size(zy_proj,3),size(zy_proj,4)];

zy_proj_2 = zeros(newsize);
zy_proj_2(1,:,:,:) = imresize3(squeeze(zy_proj(1,:,:,:)),newsize(2:end)); 
zy_proj_2(2,:,:,:) = imresize3(squeeze(zy_proj(2,:,:,:)),newsize(2:end));

xy_zy_montage = zeros(Nchan, size(xy_proj,2), size(xy_proj,3)+size(zy_proj_2,3)+5,size(xy_proj,4));

xy_zy_montage(:,:,1:size(xy_proj,3),:) = xy_proj;
xy_zy_montage(:,:,end-size(zy_proj_2,3)+1:end,:) = zy_proj_2;
%show xy/zy montage as movie
implay2chan(xy_zy_montage);
%% save xy/zy montage as tif
savepath = strcat(pipe.lab.rundir(mouse,date,run,server),'_xyzy_montage.tif');
write2chanTiff(uint16(xy_zy_montage),savepath);
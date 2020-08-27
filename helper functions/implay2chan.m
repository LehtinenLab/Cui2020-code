function rgbmov = implay2chan(mov,order)
    if nargin == 1
        order = [1,2,0];
    end
    
    if isempty(order)
        order = [1,2,0];
    end
    
    rgb_mov = zeros(size(mov,2),size(mov,3),3,size(mov,4));
    %get first channel
    rgb_mov(:,:,find(order==1),:) = rescale(squeeze(mov(1,:,:,:)));
    %get green channel
    rgb_mov(:,:,find(order==2),:) = rescale(squeeze(mov(2,:,:,:)));
    %
    implay(rgb_mov);
end


clear all

imgbk = imread('View_001\\frame_0000.jpg');

thr = 4;
minArea = 200;

baseNum = 0000;
seqLength = 794;

% baseNum = 1374;
% seqLength = 0;
% 
%imshow(imgdif)
se = strel('disk',3);
blur = [-1,0,1:-1,0,1:-1,0,1];

figure;
for i=0:seqLength
    imgfr = imread(sprintf('View_001\\frame_%.4d.jpg',baseNum+i));
    hold off
    subplot(2,2,1);imshow(imgfr);
    
    imgbkB = imfilter(imgbk, blur, 'conv');
    imgfrB = imfilter(imgfr, blur, 'conv');
    
    imgdif = (abs(double(imgbkB(:,:,1))-double(imgfrB(:,:,1)))>thr) | ...
        (abs(double(imgbkB(:,:,2))-double(imgfrB(:,:,2)))>thr) | ...
        (abs(double(imgbkB(:,:,3))-double(imgfrB(:,:,3)))>thr);
    
    imgbk = imgfr;
    
    bw = imopen(imgdif,se);
    
    subplot(2,2,2);imshow(imgbkB); 
    subplot(2,2,3);imshow(bw); drawnow
 
    [lb, num]=bwlabel(bw);
    regionProps = regionprops(lb,'area','FilledImage','Centroid');
    inds = find([regionProps.Area]>minArea);
    
    regnum = length(inds);
    
    if regnum
        for j=1:regnum
            [lin, col]= find(lb == inds(j));
            upLPoint = min([lin col]);
            dWindow  = max([lin col]) - upLPoint + 1;
           
            rectangle('Position',[fliplr(upLPoint) fliplr(dWindow)],'EdgeColor',[1 1 0],...
                'linewidth',2);
        end
    end
    drawnow
end
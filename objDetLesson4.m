

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
edgeLR = [-1,0,1:-1,0,1:-1,0,1];
edgeRL = [1,0,-1:1,0,-1:1,0,-1];

figure;
for i=0:seqLength
    imgfr = imread(sprintf('View_001\\frame_%.4d.jpg',baseNum+i));
    hold off
    subplot(2,2,1);imshow(imgfr);
    
    imgbkELR = imfilter(imgbk, edgeLR, 'conv');
    imgfrELR = imfilter(imgfr, edgeLR, 'conv');
    imgbkERL = imfilter(imgbk, edgeRL, 'conv');
    imgfrERL = imfilter(imgfr, edgeRL, 'conv');
    
    imgdif = (abs(double(imgbk(:,:,1))-double(imgfr(:,:,1)))>thr) | ...
        (abs(double(imgbk(:,:,2))-double(imgfr(:,:,2)))>thr) | ...
        (abs(double(imgbk(:,:,3))-double(imgfr(:,:,3)))>thr);
    
    imgdifLR = (abs(double(imgbkELR(:,:,1))-double(imgfrELR(:,:,1)))>thr) | ...
        (abs(double(imgbkELR(:,:,2))-double(imgfrELR(:,:,2)))>thr) | ...
        (abs(double(imgbkELR(:,:,3))-double(imgfrELR(:,:,3)))>thr);
    
    imgdifRL = (abs(double(imgbkERL(:,:,1))-double(imgfrERL(:,:,1)))>thr) | ...
        (abs(double(imgbkERL(:,:,2))-double(imgfrERL(:,:,2)))>thr) | ...
        (abs(double(imgbkERL(:,:,3))-double(imgfrERL(:,:,3)))>thr);
    
    imgbk = imgfr;
    
    bw = imopen(imgdif,se);
    bwLR = imopen(imgdifLR,se);
    bwRL = imopen(imgdifRL,se);
    
    imgFinal = bwLR + bwRL;
    
    subplot(2,2,1);imshow(imgbk);
    subplot(2,2,1);imshow(imgfr); 
    subplot(2,2,3);imshow(imgFinal);
    
    [lb, num]=bwlabel(imgFinal);
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
    
    subplot(2,2,4);imshow(bw); drawnow
 
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
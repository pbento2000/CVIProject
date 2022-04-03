

clear all

sampleXMLfile = 'PETS2009-S2L1.xml';
mlStruct = parseXML(sampleXMLfile);

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
for i=1:seqLength
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
    subplot(2,2,2);imshow(imgfr); 
    subplot(2,2,3);imshow(imgFinal);
    
    [lb, num]=bwlabel(imgFinal);
    regionProps = regionprops(lb,'area','FilledImage','Centroid');
    inds = find([regionProps.Area]>minArea);
    
    regnum = length(inds);
    genBoxes = zeros(regnum, 4);

    if regnum
        for j=1:regnum
            [lin, col]= find(lb == inds(j));
            disp([lin, col]);
            upLPoint = min([lin col]);
            dWindow  = max([lin col]) - upLPoint + 1;
            
            genBoxes(j, 1) = upLPoint(2);
            genBoxes(j, 2) = upLPoint(1);
            genBoxes(j, 3) = dWindow(2);
            genBoxes(j, 4) = dWindow(1);
            
            rectangle('Position',[fliplr(upLPoint) fliplr(dWindow)],'EdgeColor',[1 1 0],...
                'linewidth',2);
        end
    end
    
    a = size(mlStruct.Children((i+1)*2).Children(2).Children);
    a = int64((a(2)/2)-0.5);
    disp(a);
    bBoxes = zeros(a,4);
    
    for n = 1:a
        bBoxes(n,1) = str2double(mlStruct.Children((i+1)*2).Children(2).Children(n*2).Children(2).Attributes(1).Value);
        bBoxes(n,2) = str2double(mlStruct.Children((i+1)*2).Children(2).Children(n*2).Children(2).Attributes(2).Value);
        bBoxes(n,3) = str2double(mlStruct.Children((i+1)*2).Children(2).Children(n*2).Children(2).Attributes(3).Value);
        bBoxes(n,4) = str2double(mlStruct.Children((i+1)*2).Children(2).Children(n*2).Children(2).Attributes(4).Value);
    end
    
    for j=1:a
        upLPoint = [bBoxes(j, 3) - (bBoxes(j, 2)/2), bBoxes(j, 4) - (bBoxes(j, 1)/2)];
        dWindow  = [bBoxes(j, 1), bBoxes(j, 2)];
           
        bBoxes(j, 1) = upLPoint(1);
        bBoxes(j, 2) = upLPoint(2);
        bBoxes(j, 3) = dWindow(2);
        bBoxes(j, 4) = dWindow(1);
    end
    drawnow
%     subplot(2,2,4);imshow(bw); drawnow
%  
%     [lb, num]=bwlabel(bw);
%     regionProps = regionprops(lb,'area','FilledImage','Centroid');
%     inds = find([regionProps.Area]>minArea);
%     
%     regnum = length(inds);
%     
%     if regnum
%         for j=1:regnum
%             [lin, col]= find(lb == inds(j));
%             upLPoint = min([lin col]);
%             dWindow  = max([lin col]) - upLPoint + 1;
%             rectangle('Position',[fliplr(upLPoint) fliplr(dWindow)],'EdgeColor',[1 1 0],...
%                 'linewidth',2);
%         end
%     end
%     
%     drawnow
    
end
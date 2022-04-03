

clear all

sampleXMLfile = 'PETS2009-S2L1-cropped.xml';
mlStruct = parseXML(sampleXMLfile);

imgbk = imread('View_001\\frame_0000.jpg');

thr = 4;
minArea = 400;

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
           
            rectangle('Position',[upLPoint fliplr(dWindow)],'EdgeColor',[1 1 0],...
                'linewidth',2);
        end
    drawnow
end
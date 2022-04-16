clear all

sampleXMLfile = 'PETS2009-S2L1.xml';
mlStruct = parseXML(sampleXMLfile);

imgbk = imread('View_001\frame_0000.jpg');

thresholdMatrix = (0:0.05:1);
iouMatrix = zeros(1,length(thresholdMatrix));
heatmapImg = size(imgbk);
heatmapMatrix = zeros(heatmapImg(1), heatmapImg(2));

thr = 4;
minArea = 2000;

baseNum = 0000;
seqLength = 794;

step = 1;

boxId = 0;

inds_old = 0;
regionProps_old = 0;
genBoxes_old = 0;

trajectory_x = {};
trajectory_y = {};
trajectory_tail_lenght = 3;

% baseNum = 1374;
% seqLength = 0;
% 
%imshow(imgdif)
se = strel('disk',3);
edgeLR = [-1,0,1:-1,0,1:-1,0,1];
edgeRL = [1,0,-1:1,0,-1:1,0,-1];

figure;
for i=1:step:seqLength
    imgfr = imread(sprintf('View_001\\frame_%.4d.jpg',baseNum+i));
    hold off
    %subplot(1,2,1);imshow(imgfr);
    
    %Utilizar edge matrix para comparar s� as edges de cada frame
    
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
    
    bw = imclose(imgdif,se);
    bwLR = imopen(imgdifLR,se);
    bwRL = imopen(imgdifRL,se);
    
    imgFinal = bwLR + bwRL;
    %imgFinal = bw;
    
    
    subplot(1,2,1);imshow(imgfr); 
    
    %subplot(2,2,1);imshow(imgbkg);
    %subplot(2,2,2);imshow(imgFinal);;
    
    [lb, num]=bwlabel(imgFinal);
    regionProps = regionprops(lb,'area','FilledImage','Centroid');
    inds = find([regionProps.Area]>minArea);
    
    regnum = length(inds);
    genBoxes = zeros(regnum, 4);
    areaSum = 0;
    areaGTSum = 0;

    if regnum
        for j=1:regnum
            [lin, col]= find(lb == inds(j));
            pixelsRegion = [lin, col];
            centroid = regionProps(inds(j)).Centroid;
            upLPoint = min([lin col]);
            raio = sqrt((upLPoint(2) - centroid(1)).^2 + (upLPoint(1) - centroid(2)).^2);
            %disp(raio);
            dWindow  = max([lin col]) - upLPoint + 1;
            %Calcular a distancia Gaussiana
            for x = 1:length(pixelsRegion)
                %disp(pixelsRegion(x,1))
                %disp(pixelsRegion(x,2))
                heatmapMatrix(pixelsRegion(x,1),pixelsRegion(x,2)) = heatmapMatrix(pixelsRegion(x,1),pixelsRegion(x,2)) + (1-(sqrt((pixelsRegion(x,1) - centroid(2)).^2 + (pixelsRegion(x,2) - centroid(1)).^2)/raio)).^2;
            end
            %Guardar cada bounding box gerado pelo algoritmo numa matriz
            
            genBoxes(j, 1) = upLPoint(2);
            genBoxes(j, 2) = upLPoint(1);
            genBoxes(j, 3) = dWindow(2);
            genBoxes(j, 4) = dWindow(1);
            %Soma das bounding boxes calculadas
            
            areaSum = areaSum + dWindow(1) * dWindow(2);
            
            rectangle('Position',[fliplr(upLPoint) fliplr(dWindow)],'EdgeColor',[1 1 0],...
                'linewidth',2);
            
            if i > 1 && j <= length(inds_old)
                intersection_area = rectint(genBoxes_old, genBoxes);
                foundMatch = false;
                regnum_old = length(inds_old); 
                for k=1:regnum_old
                    area1 = genBoxes_old(k,3) * genBoxes_old(k,4);
                    area2 = genBoxes(j,3) * genBoxes(j,4);
                    iou = intersection_area(k,j) / (area1 + area2 - intersection_area(k,j));
                    if iou >= 0.7
                        id = genBoxes_old(k,5);
                        trajectory_x{id} = [trajectory_x{id} centroid(1)];
                        trajectory_y{id} = [trajectory_y{id} centroid(2)];
                        genBoxes(j,5) = id;
                        text(upLPoint(2),upLPoint(1),string(id),'Color','red');
                        foundMatch=true;
                        break;
                    end
                end
                if ~foundMatch
                    boxId=boxId+1;
                    trajectory_x{boxId} = centroid(1);
                    trajectory_y{boxId} = centroid(2);
                    genBoxes(j,5) = boxId;
                    text(upLPoint(2),upLPoint(1),string(boxId),'Color','red');
                end
            else
                boxId=boxId+1;
                trajectory_x{boxId} = centroid(1);
                trajectory_y{boxId} = centroid(2);
                genBoxes(j,5) = boxId;
                text(upLPoint(2),upLPoint(1),string(boxId),'Color','red');
            end
        end
    end
    
    %Desenhar trajectorias
    
    hold on;
    if boxId
        for j=1:boxId
            plot(trajectory_x{j}, trajectory_y{j}, '-cd');
            if ~mod(i,trajectory_tail_lenght)
                trajectory_x{j} = trajectory_x{j}(2:end);
                trajectory_y{j} = trajectory_y{j}(2:end);
            end
        end
    end
    hold off;
    drawnow
    
    %Atualizar informacao da frame anterior
    inds_old = inds;
    regionProps_old = regionProps;
    genBoxes_old = genBoxes;
    
    %Ler as bounding boxes do GT para este frame
    
    a = size(mlStruct.Children((i+1)*2).Children(2).Children);
    a = int64((a(2)/2)-0.5);
    %disp(a);
    bBoxes = zeros(a,4);
    
    %Guardar numa matriz
    
    for n = 1:a
        bBoxes(n,1) = str2double(mlStruct.Children((i+1)*2).Children(2).Children(n*2).Children(2).Attributes(1).Value);
        bBoxes(n,2) = str2double(mlStruct.Children((i+1)*2).Children(2).Children(n*2).Children(2).Attributes(2).Value);
        bBoxes(n,3) = str2double(mlStruct.Children((i+1)*2).Children(2).Children(n*2).Children(2).Attributes(3).Value);
        bBoxes(n,4) = str2double(mlStruct.Children((i+1)*2).Children(2).Children(n*2).Children(2).Attributes(4).Value);
        areaGTSum = areaGTSum + bBoxes(n,1) * bBoxes(n,2);
    end
    
    %Converter para o formato bounding box
    
    for j=1:a
        upLPoint = [bBoxes(j, 3) - (bBoxes(j, 2)/2), bBoxes(j, 4) - (bBoxes(j, 1)/2)];
        dWindow  = [bBoxes(j, 1), bBoxes(j, 2)];
           
        bBoxes(j, 1) = upLPoint(1);
        bBoxes(j, 2) = upLPoint(2);
        bBoxes(j, 3) = dWindow(2);
        bBoxes(j, 4) = dWindow(1);
    end
    
    area = rectint(bBoxes, genBoxes);
    sumAreaIntersection = sum(sum(area));
    
    %disp(sumAreaIntersection);
    %disp(areaSum);
    %disp(areaGTSum);
    
    %Calcular a taxa de sucesso
    
    iou = sumAreaIntersection / (areaSum + areaGTSum - sumAreaIntersection);
    
    %Guardar para diferentes thresholds de modo a poder desenhar o gr�fico
    
    for j=1:21
        if iou > thresholdMatrix(1,j)
            iouMatrix(1,j) = iouMatrix(1,j) + 1;
        end
    end
    
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
    %disp(iouMatrix);
    
    %Esta a criar um grafico para cada frame, mas pode se fazer so no fim
    subplot(1,2,2); plot(thresholdMatrix, iouMatrix,'m--o'); drawnow
end
figure;
h = heatmap(heatmapMatrix, 'Colormap', jet, 'GridVisible','off');
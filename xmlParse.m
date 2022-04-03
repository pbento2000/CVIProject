sampleXMLfile = 'PETS2009-S2L1-cropped.xml';
mlStruct = parseXML(sampleXMLfile);

nFrames = 795;

for i = 1:nFrames
    a = size(mlStruct.Children(i*2).Children(2).Children);
    a = int64((a(2)/2)-0.5);
    disp(a);
    for j = 1:a
        disp(mlStruct.Children(i*2).Children(2).Children(j*2).Children(2).Attributes(1).Value);
    end
end
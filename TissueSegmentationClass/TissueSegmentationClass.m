classdef TissueSegmentationClass < handle
    %TISSUESEGMENTATIONCLASS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        %% trial parameters 
        trialNumber;
        sliceNumber;
        numSlices; 
        needleTrackSlice; 
        needleRelativeLoc;
        sliceWidth; 
        startHeight; 
        endHeight; 
        needleHeight;
        imageList; 
        imageFileNameList;
        tifObj
        %% output data
        slices; % cell with slices
        extendedSlices; 
        masks;
        outputMask; % final output
        flipped;
        sliceInds;
        roiXcells;
        roiYcells;
        sigma;
        %% geometry
        elLength; % elevation length (mm), x, index 3
        azLength; % azimuth length (mm), y, index 2
        raLength; % range length (mm), z, index 1
        dEl;
        %% file io
        targetFolder;
        targetFile;
        currImage; 
        fileType; 
        filesInFolder;
        %% image parameters
        voxelSize; %% pixels per MM
        scaleFactor;
        imRes;
        sliceSizeInImage;
        displayExtension;
        minRange;
        maxRange;
        minAzimuth;
        maxAzimuth;
        needleLoc;
        needleLocAZRA;
        %%
        currSubImage;
    end
    
    methods
        %% Setters
        function obj = TissueSegmentationClass(raLength, azLength, elLength,voxelSize,needleLoc)
            %TISSUESEGMENTATIONCLASS Construct an instance of this class
            %   Class for tissue segmentation, contains methods and
            %   properties for use by the segmentation GUI
            obj.raLength = raLength; 
            obj.azLength = azLength; 
            obj.elLength = elLength; 
            obj.fileType = '*.tif';
            obj.voxelSize = voxelSize;
            obj.displayExtension = 10; 
            obj.needleLocAZRA = needleLoc;
            obj.sliceInds=[];
            obj.masks = {};
            obj.slices = {};
            obj.imageList = {};
            obj.imageFileNameList = {};
            obj.filesInFolder = {};
        end
        function set3DParams(obj,startHeight,endHeight, needleTrackSlice,needleRelativeLoc,needleHeight,sigma)
            % Set parameters for final mask creation
            obj.startHeight = startHeight;
            obj.endHeight = endHeight;
            obj.needleTrackSlice = needleTrackSlice;
            obj.needleHeight = needleHeight;
            obj.sigma = sigma;
            switch needleRelativeLoc
                case 'surface'
                    obj.needleRelativeLoc = 's';
                case 'middle'
                    obj.needleRelativeLoc = 'm';
                case 'bottom'     
                    obj.needleRelativeLoc = 'b';
            end
        end
        function setFolder(obj,targetFolder)
            % Set folder to pull images from and pull all .tif files
            obj.targetFolder = targetFolder; 
            % get files in folder
            obj.filesInFolder = dir(fullfile(targetFolder,'*.tif'));
            % remove temp files generated from AFS (I think) with ._ prefix, might need
            % localization for other platforms/file systems
            obj.filesInFolder=obj.filesInFolder(cellfun(@(x)~strcmp('._',x(1:2)),{obj.filesInFolder.name}));
            for currFile = 1:length(obj.filesInFolder)
                fileTarget = fullfile(obj.filesInFolder(currFile).folder,obj.filesInFolder(currFile).name);
                obj.tifObj = Tiff(fileTarget);
                obj.imRes = obj.tifObj.getTag('XResolution')/25.4;
                obj.scaleFactor = obj.voxelSize/obj.imRes; 
                obj.sliceSizeInImage = ceil([obj.imRes*obj.raLength,obj.imRes*obj.azLength]);
                obj.minRange = obj.sliceSizeInImage(1)+obj.displayExtension*obj.imRes-100;
                obj.maxRange = obj.sliceSizeInImage(1)+obj.tifObj.getTag('ImageLength')-obj.displayExtension*obj.imRes;
                obj.minAzimuth = obj.sliceSizeInImage(2)+obj.displayExtension*obj.imRes-100;
                obj.maxAzimuth = obj.sliceSizeInImage(2)+obj.tifObj.getTag('ImageWidth')-obj.displayExtension*obj.imRes;
                obj.currImage = fliplr(padarray(imread(fileTarget),obj.sliceSizeInImage));
                obj.imageList{end+1} = obj.currImage; 
                obj.imageFileNameList{end+1} = obj.filesInFolder(currFile).name; 
            end
        end
        function updateInitQuantities(obj)
            if isempty(obj.voxelSize)
                obj.voxelSize = 1; 
            end
            obj.scaleFactor = obj.voxelSize/obj.imRes; 
            obj.sliceSizeInImage = ceil([obj.imRes*obj.raLength,obj.imRes*obj.azLength]);
            obj.minRange = obj.sliceSizeInImage(1)+obj.displayExtension*obj.imRes-100;
            
            obj.minAzimuth = obj.sliceSizeInImage(2)+obj.displayExtension*obj.imRes-100;
            try
                obj.maxAzimuth = obj.sliceSizeInImage(2)+obj.tifObj.getTag('ImageWidth')-obj.displayExtension*obj.imRes;
                obj.maxRange = obj.sliceSizeInImage(1)+obj.tifObj.getTag('ImageLength')-obj.displayExtension*obj.imRes;
            catch
            end
        end
        function succ = guiSetFolder(obj)
            % Use GUI to call setfolder
            folderSelection = uigetdir('Select Target Folder'); 
            if folderSelection ~= 0
                obj.setFolder(folderSelection);
                succ = true;
            else
                succ = false;
            end
        end
        
        function succ = uiSetImage(obj)
            [obj.targetFile,filePath] = uigetfile(obj.fileType,'Select Image File',obj.targetFolder);
            obj.targetFolder = filePath;
            if filePath ~= 0
                obj.tifObj = Tiff(fullfile(filePath,obj.targetFile));
                obj.imRes = obj.tifObj.getTag('XResolution')/25.4;
                obj.scaleFactor = obj.voxelSize/obj.imRes; 
                obj.sliceSizeInImage = ceil([obj.imRes*obj.raLength,obj.imRes*obj.azLength]);
                obj.minRange = obj.sliceSizeInImage(1)+obj.displayExtension*obj.imRes;
                obj.maxRange = obj.sliceSizeInImage(1)+obj.tifObj.getTag('ImageLength')-obj.displayExtension*obj.imRes;
                obj.minAzimuth = obj.sliceSizeInImage(2)+obj.displayExtension*obj.imRes;
                obj.maxAzimuth = obj.sliceSizeInImage(2)+obj.tifObj.getTag('ImageWidth')-obj.displayExtension*obj.imRes;
                obj.currImage = fliplr(padarray(imread(fullfile(filePath,obj.targetFile)),obj.sliceSizeInImage));
                obj.imageList{end+1} = obj.currImage; 
                obj.imageFileNameList{end+1} = obj.targetFile; 
                succ = true;
            else
                succ = false; 
            end
        end
        %%
        function [subImDisplay,subIm,boxCoordsX,boxCoordsY,needlePoint] = getSubIm(obj, startVec,rotAmt,flipped)
            % Get a portion of the given RGB image
            % startVec: lower [x,y] coordinates of desired subimage
            % rotAmt: Amount to rotate in degrees
            % flipped: whether or not the image should be flipped
            if ~flipped
                rotatedImage = imrotate(obj.currImage,rotAmt);
            else
                rotatedImage = imrotate(fliplr(obj.currImage),rotAmt);
            end
            raTop = startVec(1);
            raBottom = startVec(1) + obj.sliceSizeInImage(1); 
            azTop = startVec(2);
            azBottom = startVec(2) + obj.sliceSizeInImage(2); 
            raVec = raTop:raBottom;
            azVec = azTop:azBottom;
            raVecExt = (raTop-ceil(obj.imRes*obj.displayExtension)):(raBottom+ceil(obj.imRes*obj.displayExtension));
            azVecExt = (azTop-ceil(obj.imRes*obj.displayExtension)):(azBottom+ceil(obj.imRes*obj.displayExtension));
            subIm = rotatedImage(raVec,azVec,:);
            subImDisplay = rotatedImage(raVecExt,azVecExt,:);
            boxRa1 = ceil(obj.imRes*obj.displayExtension);
            boxRa2 = ceil(obj.imRes*obj.displayExtension)+obj.sliceSizeInImage(1);
            boxAz1 = ceil(obj.imRes*obj.displayExtension);
            boxAz2 = ceil(obj.imRes*obj.displayExtension)+obj.sliceSizeInImage(2);
            boxCoordsX = [boxRa1,boxRa2,boxRa2,boxRa1,boxRa1];
            boxCoordsY = [boxAz1,boxAz1,boxAz2,boxAz2,boxAz1];
            needlePoint = [boxRa1,boxAz1]+obj.needleLocAZRA*obj.imRes;
        end
        function [subImDisplay,subIm,boxCoordsX,boxCoordsY,needlePoint] = getSubImInd(obj, ind, startVec,rotAmt,flipped)
            if ~flipped
                rotatedImage = imrotate(obj.imageList{ind},rotAmt);
            else
                rotatedImage = imrotate(fliplr(obj.imageList{ind}),rotAmt);
            end
            raTop = startVec(1);
            raBottom = startVec(1) + obj.sliceSizeInImage(1); 
            azTop = startVec(2);
            azBottom = startVec(2) + obj.sliceSizeInImage(2); 
            raVec = raTop:raBottom;
            azVec = azTop:azBottom;
            raVecExt = (raTop-ceil(obj.imRes*obj.displayExtension)):(raBottom+ceil(obj.imRes*obj.displayExtension));
            azVecExt = (azTop-ceil(obj.imRes*obj.displayExtension)):(azBottom+ceil(obj.imRes*obj.displayExtension));
            subIm = rotatedImage(raVec,azVec,:);
            subImDisplay = rotatedImage(raVecExt,azVecExt,:);
            boxRa1 = ceil(obj.imRes*obj.displayExtension);
            boxRa2 = ceil(obj.imRes*obj.displayExtension)+obj.sliceSizeInImage(1);
            boxAz1 = ceil(obj.imRes*obj.displayExtension);
            boxAz2 = ceil(obj.imRes*obj.displayExtension)+obj.sliceSizeInImage(2);
            boxCoordsX = [boxRa1,boxRa2,boxRa2,boxRa1,boxRa1];
            boxCoordsY = [boxAz1,boxAz1,boxAz2,boxAz2,boxAz1];
            needlePoint = [boxRa1,boxAz1]+obj.needleLocAZRA*obj.imRes;
        end
        function funcOut = getSubImClosure(obj,fileName)
            ind = cellfun(@(x)strcmp(x,fileName), obj.imageFileNameList);
            funcOut = @(obj, startVec,rotAmt,flipped) obj.getSubImInd(ind, startVec,rotAmt,flipped);
        end
        function addSlice(obj, ind, dat,datPadded,flipped)
            obj.slices{end+1} = dat; 
            obj.extendedSlices{end+1}=datPadded;
            obj.sliceInds(end+1) = ind;
            obj.flipped(end+1) = flipped; 
        end
        function modifySlice(obj, ind, dat,datPadded,flipped)
            slInd = find(obj.sliceInds == ind);
            obj.slices{slInd} = dat; 
            obj.extendedSlices{slInd}=datPadded;
            obj.sliceInds(slInd) = ind;
            obj.flipped(slInd) = flipped; 
        end
        function sliceExists = sliceExists(obj, ind)
             sliceExists = ~isempty((find(obj.sliceInds==ind)));
        end
        function maskExists = maskExists(obj, ind)
            slInd = find(obj.sliceInds==ind); 
            try
                if length(slInd)~=0
                    maskExists = true;
                else
                    maskExists = false;
                end
            catch
                maskExists = false;
            end
        end
        function outIm = getSlice(obj, ind) 
            outIm = obj.slices{obj.sliceInds == ind};
        end
        function outIm = getSlicePadded(obj, ind)
            outIm = obj.extendedSlices{obj.sliceInds == ind};
        end
        function outIm = getMask(obj, ind)
            outIm = obj.masks{obj.sliceInds == ind};
        end
        function outIm = getMaskPadded(obj, ind)
            outIm = padarray(obj.masks{obj.sliceInds == ind},[ceil(obj.imRes*obj.displayExtension),ceil(obj.imRes*obj.displayExtension)]);
        end
        function [roiX,roiY] = getMaskPoly(obj, ind)
            slInd = obj.sliceInds == ind;
            offset = ceil(obj.imRes*obj.displayExtension); 
            roiX = obj.roiXcells{slInd}+offset;
            roiY = obj.roiYcells{slInd}+offset;
        end
        function [roiX,roiY] = getMaskPolyUnpadded(obj, ind)
            roiX = obj.roiXcells{obj.sliceInds == ind};
            roiY = obj.roiYcells{obj.sliceInds == ind};
        end
        function padAmt = getPadding(obj)
            padAmt = ceil(obj.imRes*obj.displayExtension);
        end
        function createMask(obj, ind)
            slInd = find(obj.sliceInds == ind);
            fig = figure; 
            imshow(obj.slices{slInd}); axis image;
            set(fig, 'MenuBar', 'none');
            set(fig, 'ToolBar', 'none');
            [roi,roiX,roiY] = roipoly();
            close(fig);
            obj.masks{slInd} = roi;
            obj.roiXcells{slInd} = roiX;
            obj.roiYcells{slInd} = roiY;
        end
        function createMaskFromPoly(obj, ind,polyIn, maskIn)
            slInd = find(obj.sliceInds == ind);
            roiX = ceil(polyIn(:,1)');
            roiY = ceil(polyIn(:,2)');
            roiX = [roiX, roiX(1)];
            roiY = [roiY, roiY(1)];
            obj.masks{slInd} = maskIn;
            obj.roiXcells{slInd} = roiX;
            obj.roiYcells{slInd} = roiY;
        end
        function outputMask = create3DVolumeOld(obj,interpMethod,adjCenter)
            [sortedInds,I] = sort(obj.sliceInds);
            sortedMasks = obj.masks(I);
            sortedFlipped = obj.flipped(I);
            maskConcat = (cat(3,sortedMasks{:}));
            regionPixels = ceil(obj.imRes*(obj.endHeight-obj.startHeight));
            obj.dEl = regionPixels/length(obj.masks);
            elVec = zeros(1,length(sortedMasks));
            for currMask = 1:length(sortedMasks)
                locStart = floor(obj.dEl*(currMask-1))+1;
                locEnd = floor(obj.dEl*currMask);
                if ~sortedFlipped(currMask)
                    elVec(currMask) = locEnd;
                else
                    elVec(currMask) = locStart;
                end
            end
            raVec = 1:size(sortedMasks{1},1);
            azVec = 1:size(sortedMasks{1},2);
            [raMat,azMat,elMat] = meshgrid(raVec,azVec,elVec);
            [raMatT,azMatT,elMatT] = meshgrid(raVec,azVec,1:regionPixels);
            if ~strcmp(interpMethod,'pchip')
                maskConcat = double(maskConcat);
                outputMask = ((interp3(raMat,azMat,elMat,maskConcat,raMatT,azMatT,elMatT,interpMethod)));
                outputMask(isnan(outputMask)) = 0;
            else
                outputMask = double(interpmask(1:size(maskConcat,3), (maskConcat), linspace(1,size(maskConcat,3),regionPixels),'pchip'));
            end
            %% needle track
            obj.needleLoc = elVec(find(sortedInds==obj.needleTrackSlice));
            if ~sortedFlipped(find(sortedInds==obj.needleTrackSlice))
                if obj.needleRelativeLoc == 'b'
                    obj.needleLoc = round(obj.needleLoc - obj.dEl);
                elseif obj.needleRelativeLoc == 'm'
                    obj.needleLoc = round(obj.needleLoc - obj.dEl/2);
                end
            else
                if obj.needleRelativeLoc == 'b'
                    obj.needleLoc = round(obj.needleLoc + obj.dEl);
                elseif obj.needleRelativeLoc == 'm'
                    obj.needleLoc = round(obj.needleLoc + obj.dEl/2);
                end
            end
            finalNeedleLoc = round(obj.needleHeight*obj.imRes); 
            bottomPadAmount = finalNeedleLoc - obj.needleLoc;
            outputMask = padarray(outputMask,[0,0,bottomPadAmount],'pre');
            currHeight = size(outputMask,3)*obj.scaleFactor;
            topPadAmount = ceil((obj.elLength-currHeight)/obj.scaleFactor);
            outputMask = padarray(outputMask,[0,0,topPadAmount],'post');
            
            if adjCenter
                [raInd,azInd,~] = ind2sub(size(outputMask),find(outputMask == 1));
                outMomentRa = mean(raInd);
                outMomentAz = mean(azInd);
                outputMask = noncircshift(outputMask, [floor((obj.needleLocAZRA(1)*obj.imRes)-floor(outMomentRa)),floor((obj.needleLocAZRA(2)*obj.imRes)-floor(outMomentAz)),0]);
            end
            outputMask = logical(round(imgaussfilt3(imresize3(outputMask,[obj.raLength,obj.azLength,obj.elLength],'nearest'),obj.sigma)));
            obj.outputMask = outputMask; 
            figure,isosurface(outputMask);
        end
        function [rangeCentroid, azimuthCentroid, elevationCentroid,volume] = getVolumeStats(obj)
            [raInd,azInd,elInd] = ind2sub(size(obj.outputMask),find(obj.outputMask == 1));
            rangeCentroid = mean(raInd);
            azimuthCentroid = mean(azInd);
            elevationCentroid = mean(elInd);
            % temp fix
            if isempty(obj.voxelSize)
                obj.voxelSize = 1; 
            end
            volume = .001*obj.voxelSize^3 * sum(obj.outputMask == 1,'all');
        end
        function deleteMask(obj,ind)
            slInd = find(obj.sliceInds == ind);
            obj.slices{slInd} = []; 
            obj.sliceInds(slInd) = [];
            obj.flipped(slInd) = [];
            obj.masks(slInd) = []; 
        end
        function [outImAz,outImRa, needleLoc,raVec,azVec,elVec] = geometryPreview(obj)
            [sortedInds,I] = sort(obj.sliceInds);
            sortedMasks = obj.masks(I);
            sortedFlipped = obj.flipped(I);
            maskConcat = double(cat(3,sortedMasks{:}));
            regionPixels = ceil(obj.imRes*(obj.endHeight-obj.startHeight));
            obj.dEl = regionPixels/length(obj.masks);
            elVec = zeros(1,length(sortedMasks));
            for currMask = 1:length(sortedMasks)
                locStart = floor(obj.dEl*(currMask-1))+1;
                locEnd = floor(obj.dEl*currMask);
                if ~sortedFlipped(currMask)
                    elVec(currMask) = locEnd;
                else
                    elVec(currMask) = locStart;
                end
            end
            raVec = 1:size(sortedMasks{1},1);
            azVec = 1:size(sortedMasks{1},2);
            obj.needleLoc = elVec(sortedInds==obj.needleTrackSlice);
            if ~sortedFlipped(find(sortedInds==obj.needleTrackSlice))
                if obj.needleRelativeLoc == 'b'
                    obj.needleLoc = round(obj.needleLoc - obj.dEl);
                elseif obj.needleRelativeLoc == 'm'
                    obj.needleLoc = round(obj.needleLoc - obj.dEl/2);
                end
            else
                if obj.needleRelativeLoc == 'b'
                    obj.needleLoc = round(obj.needleLoc + obj.dEl);
                elseif obj.needleRelativeLoc == 'm'
                    obj.needleLoc = round(obj.needleLoc + obj.dEl/2);
                end
            end
            finalNeedleLoc = round(obj.needleHeight*obj.imRes); 
            bottomPadAmount = finalNeedleLoc - obj.needleLoc;
            totalHeight = ceil((obj.elLength*obj.imRes));
            outImAz = zeros(totalHeight,length(azVec));
            outImRa = zeros(totalHeight,length(raVec));
            elVecIm = elVec+bottomPadAmount;
            for n = 1:3
                elVecIm = [elVecIm,elVec+bottomPadAmount+n];
            end
            outImAz(elVecIm,:) = 1; 
            outImRa(elVecIm,:) = 1; 
            needleLoc = [finalNeedleLoc*obj.scaleFactor,ceil(obj.needleLocAZRA(1))];
            azVec = azVec * obj.scaleFactor;
            elVec = (1:totalHeight)*obj.scaleFactor;
            raVec = raVec * obj.scaleFactor; 
            obj.dEl = obj.dEl*obj.scaleFactor; 
        end
    end
    
end


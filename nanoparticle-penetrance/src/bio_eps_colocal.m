%experiments = [4 5 6 7 8 14 15 16 17 18 19 20];
experiments = [4 5 6 7 8 15 17 18 19 20];
%experiments = [17];
treatment = true;
maxStackOut = zeros(6,length(experiments)*2);
maxConcOut = zeros(6,length(experiments)*2);
%EPSOut = zeros(6,length(experiments)*2);
bioOut = zeros(6,length(experiments)*2);
colocFracOut = zeros(6,length(experiments)*2);
tic;
% apparently can't parfor
for k=1:length(experiments)
    expSet = experiments(k);
    for expRep=1:2
        outputIdx = k + (k-1) + (expRep-1);
        dirStr1 = sprintf('Pilot %d',expSet);
        if treatment == true
            dirStr2 = sprintf('Treatment %d',expRep);
            dirStr3 = sprintf('P%d_D%d',expSet, expRep);
        else
            dirStr2 = sprintf('Control %d',expRep);
            dirStr3 = sprintf('P%d_ND%d',expSet, expRep);
        end
        
        %dirStr = append("./",dirStr1,"/","2 channel .mat + BAIT output/",dirStr2,"/",dirStr3,"_GR_UE.mat");
        dirStr = append("./",dirStr1,"/","2 channel .mat + BAIT output/",dirStr2,"/",dirStr3,"_RB_UE.mat");
        %data = struct2cell(load("./Pilot 7/2 channel .mat + BAIT output/Treatment 1/P7_D1_GB.mat"));
        data = struct2cell(load(dirStr));
        data = data{1};
        R = data(:,:,:,:,1);
        B = data(:,:,:,:,2);
        dims = size(R);
        height = dims(3);
        %colormap copper
        %hold on
        out = zeros(height,6);
        maxStackVec = zeros(6,1);
        maxConcVec = zeros(6,1);
        densityVec = zeros(6,1);
        %EPSVec = zeros(6,1);
        bioVec = zeros(6,1);
        parfor t=1:6
            M = R(:,:,:,t);
            %bioMat = R(:,:,:,t);
            EPSMat = B(:,:,:,t);
            imageDim = size(M);
            imageVoxels = prod(imageDim);
            BW = M > 14;
            EPSBW = EPSMat > 14;
            %bioBW = bwareaopen(bioBW,7);
            %bioMat = bioMat.*uint8(bioBW);
            EPSBySlice = squeeze(sum(EPSMat,[1 2]));
            dx = diff(EPSBySlice);
            maxDx = max(dx);
            idxOfMax = find(dx==maxDx);
            
            %EPSBW = EPSMat > 14;
            overlap = BW & EPSBW;
            coloc = EPSMat.*uint8(overlap);
            coloc(:,:,1:idxOfMax) = 0;
            overlap(:,:,1:idxOfMax) = 0;
            EPSBW(:,:,1:idxOfMax) = 0;
            BW(:,:,1:idxOfMax) = 0;
            M(:,:,1:idxOfMax) = 0;
            %BW = bwareaopen(BW,7);
            %EPSBW = bwareaopen(EPSBW,7);

            totalBio = sum(sum(sum(M)));
            totalcoloc = sum(sum(sum(coloc)));
            bv = zeros(imageDim(3),1);
            densityVec(t) = totalcoloc/imageVoxels;
            for i=1:imageDim(3)
                bv(i) = sum(sum(coloc(:,:,i)))/totalcoloc;
            end
            out(:,t) = bv;
            maxConc = max(bv);
            maxStack = find(bv==maxConc);
            colocFracOfBio = sum(coloc,[1 2 3])/totalBio;
            maxConcVec(t) = maxConc;
            maxStackVec(t) = (height - maxStack(1))/height;
            bioVec(t) = totalBio;
            colocFracVec(t) = colocFracOfBio;
            %EPSVec(t) = totalEPS;
            %plot(bv)
        end
        maxStackOut(:,outputIdx) = maxStackVec;
        maxConcOut(:,outputIdx) = maxConcVec;
        bioOut(:,outputIdx) = bioVec;
        colocFracOut(:,outputIdx) = colocFracVec;
        %EPSOut(:,outputIdx) = EPSVec;
    end
end
toc;
%EPSPctDiff = (EPSOut(6,:)-EPSOut(1,:))./EPSOut(1,:);
EPSPctDiff = (bioOut(6,:)-bioOut(1,:))./bioOut(1,:);
avgMaxStack = median(maxStackOut,2);
lower = quantile(maxStackOut,0.1,2);
upper = quantile(maxStackOut,0.9,2);
%plot(1:6,avgMaxStack);
%hold on
%plot(upper);
%plot(lower);
%legend('median','upper','lower');
%set(gca, 'ColorOrder', copper(6))
%bar3(out)
%xlabel('Time Point')
%ylabel('Z-Stack')
%zlabel('Proportion of Total Biovolume')
%legend('t=10','t=20','t=30','t=40','t=50','t=60','Location','northeast')
%hold off
%axis square;
%set(gca, 'Ydir', 'reverse');
%figure();
%plot(densityVec)
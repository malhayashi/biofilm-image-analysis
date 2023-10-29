function [movementVec, objectVec, matName] = calculate_object_mvmt(input_struct,idx)
    S = input_struct;
    tic;
    fields = fieldnames(S);
    M = S.(fields{idx});
    matName = fields{idx};
    BW = M > 14;
    %BW = threshold_1channel_binary_time(M);
    % Get centroids for a sequential pair of frames
    T = size(BW,4);
    movementVec = zeros(T-1,1);
    objectVec = zeros(T-1,1);
    parfor t=1:T-1
        frame1 = BW(:,:,:,t);
        filteredFrame1 = bwareaopen(frame1,6,26);
        CC1 = bwconncomp(filteredFrame1);
        Centroids1 = regionprops(CC1,'Centroid');
        
        frame2 = BW(:,:,:,t+1);
        filteredFrame2 = bwareaopen(frame2,6,26);
        CC2 = bwconncomp(filteredFrame2);
        Centroids2 = regionprops(CC2,'Centroid');
        
        minDists = zeros(length(Centroids1),1);
        for i=1:length(Centroids1)
            dists = zeros(length(Centroids2),1);
            for j=1:length(Centroids2)
                dist = norm(Centroids1(i).Centroid-Centroids2(j).Centroid);
                dists(j) = dist;
            end
            minDists(i) = min(dists);
        end
        movement=sum(minDists);
        movementVec(t) = movement;
        objectVec(t) = CC2.NumObjects
    end
    toc
end

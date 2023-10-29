function [biovolume,matName] = calculate_biovolume(input_struct,idx)
    S = input_struct;
    tic;
    fields = fieldnames(S);
    M = S.(fields{idx});
    matName = fields{idx};
    BW = M > 14;
    %BW = threshold_1channel_binary_time(M);
    % Get centroids for a sequential pair of frames
    T = size(BW,4);
    biovolume = zeros(T,1);
    parfor t=1:T
        frame = BW(:,:,:,t);
        filteredFrame = bwareaopen(frame,6,26);
        biovolume(t) = sum(sum(sum(filteredFrame)));
    end
    toc
end
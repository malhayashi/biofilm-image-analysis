function [surfaceArea, matName] = calculate_surface_area(input_struct,idx)
    S = input_struct;
    tic;
    fields = fieldnames(S);
    M = S.(fields{idx});
    matName = fields{idx};
    BW = M > 14;
    %BW = threshold_1channel_binary_time(M);
    T = size(BW,4);
    surfaceArea = zeros(T,1);
    parfor t=1:T
        frame = BW(:,:,:,t);
        filteredFrame = bwareaopen(frame,6,26);
        CC = bwconncomp(filteredFrame);
        frameArea = table2array(regionprops3(CC, 'SurfaceArea'));
        surfaceArea(t) = sum(frameArea);
    end
    toc
end
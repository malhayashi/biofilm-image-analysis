function out_mat = threshold_1channel_binary_time(in_mat)
    matrix_size=size(in_mat);
    times = matrix_size(4);
    
    out_mat = zeros(matrix_size);
    parfor t = 1:times
        M = in_mat(:,:,:,t)
        threshold = threshold_1channel(M);
        tempMat = M > threshold;
        out_mat(:,:,:,t) = uint8(tempMat);
    end
end

function out_mat = threshold_1channel_time(in_mat)
    matrix_size=size(in_mat);
    times = matrix_size(4);
    
    out_mat = zeros(matrix_size);
    parfor t = 1:times
        threshold = threshold_1channel(in_mat,t);
        tempMat = in_mat(:,:,:,t) > threshold;
        tempMat = uint8(tempMat);
        out_mat(:,:,:,t) = in_mat(:,:,:,t).*tempMat;
    end
end

    
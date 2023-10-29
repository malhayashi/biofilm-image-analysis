function threshold_value = threshold_1channel(M)
    % BEM Thresholding
    
    vcols = zeros(1,256);
    
    for i = 0:255
        BW = M > i;       % Use a boolean filter directly
        vcols(i+1) = sum(BW(:));  % sum(A(:)) collapses all dimensions
    end
    
    threshold = 1:1:256;
    
    [fitObj, gof] = fit(transpose(threshold), transpose(vcols(1,:)), 'power2');
    
    range=(1:256);
    
    [fx] = differentiate(fitObj, range);
    
    percent_change=[transpose(1:255)];
    
    percent_change = [percent_change abs(diff(fx)./fx(1:end-1,:))];
    
    thresholds = [];
    
    for i = 1:255
        if percent_change(i,2) < .10
            thresholds = [thresholds percent_change(i,:)];
            break
        end
    end
    
    threshold_value = thresholds(1,1);
end
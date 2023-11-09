function Feature = ECGFeature(AllData,fs)
    [s,f,~] = fsst(AllData,fs,kaiser(128));
    f_indices = (f > 0.5) & (f < 40);
    fea1 = real(s(f_indices,:))';
    fea2 = imag(s(f_indices,:))';
    Feature = num2cell([fea1,fea2],2);
    XV = [Feature{:}];
    mu = mean(XV,2);
    sg = std(XV,[],2);
    Feature = cellfun(@(x)(x-mu)./sg,Feature,'UniformOutput',false);
end
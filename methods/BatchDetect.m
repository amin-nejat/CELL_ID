function BatchDetect(file,mp_params)
    [data, info, prefs, ~, ~] = NeuroPALImage.open([file, '.czi']);

    data_RGBW = double(data(:,:,:,prefs.RGBW));
    data_zscored_raw = Preprocess.zscore_frame(double(data_RGBW));

    [~, mask] = Preprocess.filter_gut_lysosomes(data_zscored_raw);
    mask = repmat(mask,1,1,1,length(prefs.RGBW));
    data_zscored = data_zscored_raw; data_zscored(mask) = 0;

    filter = AutoDetect.get_filter(mp_params.hnsz, mp_params.hnsz, 0);
    sp = AutoDetect.instance().detect(data_zscored, filter, mp_params.k, mp_params.min_eig_thresh, info.scale', mp_params.exclusion_radius);

    save_for_parfor([file, '_sp', '.mat'], sp, mp_params);
end

function save_for_parfor(fname,sp,mp_params)
    save(fname, 'sp', 'mp_params');
end
          

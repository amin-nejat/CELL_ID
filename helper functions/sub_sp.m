function subsp = sub_sp(sp, subset)
% Calculating a subset of superpixels.
%
% Amin Nejat

    subsp = sp;
    
    fields = fieldnames(subsp);
    
    for t = 1: length(sp)
        for field_index = 1: size(fields)
            try
                subsp(t).(fields{field_index}) = sp(t).(fields{field_index})(subset,:,:,:,:);
            catch
            end
        end
    end
    
end
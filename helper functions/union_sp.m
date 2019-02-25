function unionsp = union_sp(sp1, sp2)

% Calculating a subset of superpixels.
%
% Amin Nejat

    unionsp = sp1;
    
    fields = fieldnames(unionsp);
    
    for t = 1: length(sp1)
        for field_index = 1: size(fields)
            try
                unionsp(t).(fields{field_index})(end+1:end+size(sp2(t).(fields{field_index}),1),:,:,:,:) = sp2(t).(fields{field_index})(:,:,:,:,:);
            catch
            end
        end
    end
end
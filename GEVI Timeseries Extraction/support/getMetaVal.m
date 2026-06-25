function val = getMetaVal(metadata, field)

ind1 = strfind(metadata, field);

if (~isempty(ind1))
    ind2 = strfind(metadata(ind1:end), 10);
    ind2 = ind2(1);
    line = metadata(ind1:ind1+ind2-2);
    
    % Check if the line contains an array after '='
    eq_ind = strfind(line, '=');
    if ~isempty(eq_ind)
        bracket_ind = strfind(line(eq_ind:end), '[');
        if ~isempty(bracket_ind)
            val = str2num(line(bracket_ind(1)+eq_ind-1:end));
        else
            val = sscanf(line, [field ' = %f']);
        end
    else
        val = sscanf(line, [field ' = %f']);
    end
else
    disp('Error - Invalid metadata field');
    val = -1;
end

function [prot,rstraj] = read_twix_hdr(fid)

% function to read raw data header information from siemens MRI scanners 
% (currently VB and VD software versions are supported and tested).
%
% Author: Philipp Ehses MPI Tuebingen, Mar/11/2014
% email: philipp.ehses@dzne.de
      
    isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;
    
    nbuffers = fread(fid, 1,'uint32');
    
    prot = [];
    for b=1:nbuffers
        %now read string up to null termination     
        bufname = fread(fid, 10, 'uint8=>char').';
        if isOctave
              bufname(~isascii (bufname)) = [];   % (ND) removes non-ascii characters to fix errors in Octave
        end
        bufname = regexp(bufname, '^\w*', 'match');
        bufname = bufname{1};
        fseek(fid, numel(bufname)-9, 'cof');        
        buflen         = fread(fid, 1,'uint32');
        buffer         = fread(fid, buflen, 'uint8=>char').';
        buffer         = regexprep(buffer,'\n\s*\n',''); % delete empty lines
        prot.(bufname) = parse_buffer(buffer);
    end
    
    if nargout>1
        rstraj = [];
        if isfield(prot.Meas,'alRegridMode') && prot.Meas.alRegridMode(1)>1
            ncol      = prot.Meas.alRegridDestSamples(1);
            dwelltime = prot.Meas.aflRegridADCDuration(1)/ncol;
            gr_adc    = zeros(1,ncol,'single');
%             start     = prot.Meas.alRegridRampupTime(1) - (prot.Meas.aflRegridADCDuration(1)-prot.Meas.alRegridFlattopTime(1))/2;
            start     = prot.Meas.alRegridDelaySamplesTime(1);
            time_adc  = start + dwelltime * (0.5:ncol);            
            ixUp      = time_adc <= prot.Meas.alRegridRampupTime(1);
            ixFlat    = (time_adc <= prot.Meas.alRegridRampupTime(1)+prot.Meas.alRegridFlattopTime(1)) & ~ixUp;
            ixDn      = ~ixUp & ~ixFlat;
            gr_adc(ixFlat) = 1;
            if prot.Meas.alRegridMode(1) == 2  
                % trapezoidal gradient
                gr_adc(ixUp)   = time_adc(ixUp)/prot.Meas.alRegridRampupTime(1);
                gr_adc(ixDn)   = 1 - (time_adc(ixDn)-prot.Meas.alRegridRampupTime(1)-prot.Meas.alRegridFlattopTime(1))/prot.Meas.alRegridRampdownTime(1);
            elseif prot.Meas.alRegridMode(1) == 4  
                % sinusoidal gradient
                gr_adc(ixUp)   = sin(pi/2*time_adc(ixUp)/prot.Meas.alRegridRampupTime(1));
                gr_adc(ixDn)   = sin(pi/2*(1+(time_adc(ixDn)-prot.Meas.alRegridRampupTime(1)-prot.Meas.alRegridFlattopTime(1))/prot.Meas.alRegridRampdownTime(1)));
            else
                warning('regridding mode unknown');
                return;
            end
            % make sure that gr_adc is always positive (rstraj needs to be
            % strictly monotonic):
            gr_adc = max(gr_adc, 1e-4);
            rstraj = (cumtrapz(gr_adc(:)) - ncol/2)/sum(gr_adc(:));
            rstraj = rstraj - mean(rstraj(ncol/2:ncol/2+1));
            % scale rstraj by kmax (only works if all slices have same FoV!!!)
            kmax = prot.MeasYaps.sKSpace.lBaseResolution/...
                prot.MeasYaps.sSliceArray.asSlice{1}.dReadoutFOV;
            rstraj = kmax * rstraj;
        end
    end
    
end        


function prot = parse_buffer(buffer)
    [ascconv, xprot] = regexp(buffer,'### ASCCONV BEGIN[^\n]*\n(.*)\s### ASCCONV END ###','tokens','split');

    if ~isempty(ascconv)
        ascconv = [ascconv{:}{:}];
        prot = parse_ascconv(ascconv);
    else
        prot = struct();
    end

    if ~isempty(xprot)
        xprot = strcat(xprot{:}); % bug fix by Qiuting Wen (added strcat)
        xprot = parse_xprot(xprot);
        if isstruct(xprot)
            name   = cat(1,fieldnames(prot),fieldnames(xprot));
            val    = cat(1,struct2cell(prot),struct2cell(xprot));
            [~,ix] = unique(name);
            prot   = cell2struct(val(ix),name(ix));
        end
    end
end


function xprot = parse_xprot(buffer)
    xprot = [];
    tokens = regexp(buffer, '<Param(?:Bool|Long|String)\."(\w+)">\s*{([^}]*)','tokens');
    tokens = [tokens, regexp(buffer, '<ParamDouble\."(\w+)">\s*{\s*(<Precision>\s*[0-9]*)?\s*([^}]*)','tokens')];
    for m=1:numel(tokens)
        name         = char(tokens{m}(1));
        % field name has to start with letter
        if (~isletter(name(1)))
            name = strcat('x', name);
        end

        value = char(strtrim(regexprep(tokens{m}(end), '("*)|( *<\w*> *[^\n]*)', '')));
        value = regexprep(value, '\s*', ' ');

        try %#ok<TRYNC>
            value = eval(['[' value ']']);  % inlined str2num()
        end

        xprot.(name) = value;
    end
    %% ParamArrays 
    ind = findstr(buffer,'ParamArray');
    for i = 1:length(ind) 
        if i < length(ind)
            next = i+1;
            while (ind(next) < ind(i)+5000 & next < length(ind)), next = next+1; end;
            stubArr = buffer((ind(i)-1):min(length(buffer),ind(next)-1));
        else
            stubArr = buffer((ind(i)-1):length(buffer));
        end;
        namet = regexp(stubArr, '<ParamArray\."(\w+)">','tokens');
        if (~isempty(namet))
            workarr = stubArr;
            name = safe_fieldname(namet{1}{1});
            stubArr = extract_brace_string(stubArr);
            [tagStr2, tagType, stubArr] = find_next_tag(stubArr);
            if (~strcmpi(tagStr2,'Default'))
                stubArr = extract_brace_string(workarr);
            end
            tmp = [];level = 0;
            tmp = parse_loop(tmp, stubArr, namet, level);

            if (~isempty(name) & ~isfield(xprot,name))
                xprot.(name) = tmp.x;
            end;
        end;
    end;
end


function mrprot = parse_ascconv(buffer)  
    mrprot = [];    
    % [mv] was: vararray = regexp(buffer,'(?<name>\S*)\s*=\s(?<value>\S*)','names');
    vararray = regexp(buffer,'(?<name>\S*)\s*=\s*(?<value>\S*)','names');
    
    isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;
    
    for var=vararray

        try
            value = eval(['[' var.value ']']);  % inlined str2num()
        catch
            value = var.value;
        end
        
        % now split array name and index (if present)
        v = regexp(var.name,'(?<name>\w*)\[(?<ix>[0-9]*)\]|(?<name>\w*)','names');

        cnt = 0;
        tmp = cell(2, numel(v));

        breaked = false;
        for k=1:numel(v)
            if isOctave
                %vk = v{k};  % (ND) Octave error: struct cannot be indexed with {
                vk = v(k);
                if iscell(vk.name)
                    % lazy fix that throws some info away
                    vk.name = vk.name{1};
                    vk.ix   = vk.ix{1};
                end
            else
                vk = v(k);
            end
            if ~isletter(vk.name(1))
                breaked = true;
                break;
            end
            cnt = cnt+1;
            tmp{1,cnt} = '.';
            tmp{2,cnt} = vk.name;

            if ~isempty(vk.ix)
                cnt = cnt+1;
                tmp{1,cnt} = '{}';
                tmp{2,cnt}{1} = 1 + str2double(vk.ix);
            end
        end
        if ~breaked && ~isempty(tmp)
            S = substruct(tmp{:});
            mrprot = subsasgn(mrprot,S,value);
        end
    end 
end

function stvar = safe_fieldname(tagStr)
    % this function checks potential fieldnames and makes sure they are valid
    % for MATLAB syntax, e.g. 2DInterpolation -> x2DInterpolation (must begin
    % with a letter)
    
    tagStr = strtrim(tagStr);
    
    if (isletter(tagStr(1)))
        stvar = tagStr;
    else
        stvar = strcat('x', tagStr);
    end
    
    if strfind(stvar, ';'), stvar = strrep(stvar, ';', '_'); end
    if strfind(stvar, '@'), stvar = strrep(stvar, '@', '_'); end % VD13
    if strfind(stvar, '-'), stvar = strrep(stvar, '-', '_'); end
end

function stvar = extract_brace_string(text)
    % from parse_xprot.m from E. Auerbach, CMRR, 2013
    % extracts string from within curly braces, including nested braces

    tstart = strfind(text,'{');
    tend = strfind(text,'}');
    
    stvar = [];
    if (~isempty(tstart) & ~isempty(tend))
        [brackind,ind] = sort([tstart,tend]);
        pos = [ones(1,length(tstart)),-1*ones(1,length(tend))];
        endind = find(cumsum(pos(ind))==0,1);
        if (~isempty(endind))
            endpos = brackind(endind);
            stvar = text(tstart(1)+1:endpos-1);
        end;
    end;
end

function [tagStr, tagType, remStr] = find_next_tag(inStr)
    % from parse_xprot.m from E. Auerbach, CMRR, 2013
    % returns <tag> name and the remainder of the string following the tag.
    % for e.g. <Tag>, returns tagStr='Tag', tagType=''
    % for e.g. <ParamLong."Tag">, returns tagStr='Tag', tagType='ParamLong'
    % if no tag is found, returns null strings
    
    tagStr = [];
    tagType = [];
    remStr = [];
    
    startPos = strfind(inStr,'<'); % look for start of tag
    if (startPos)
        endPos = strfind(inStr,'>'); % look for end of tag
        if (endPos)
            % found complete tag
            fullTag = inStr(startPos+1:endPos-1);
            
            % now check for name/type
            dotPos = strfind(fullTag,'."');
            if (dotPos)
                tagStr = getQuotString(fullTag(dotPos+1:end));
                tagType = fullTag(1:dotPos-1);
            else
                tagStr = fullTag;
            end
            
            % return remainder
            remStr = inStr(endPos+1:end);
        end
    end
end

function mrprot = parse_loop(mrprot, workarr, tagName, level)
    % from parse_xprot.m from E. Auerbach, CMRR, 2013
    % this function is called recursively to parse the xprotocol
    
    [tagStr, tagType, workarr] = find_next_tag(workarr);
    while (~isempty(workarr))
        % for parammap, add the name as another level of the current array name
        % and spawn off another copy of this function to deal with them
        if (strncmpi(tagType, 'ParamMap', 8))
            if (~isempty(tagStr))
                level = level + 1;
                tagName{level} = safe_fieldname(tagStr);
            else
                level = level + 1;
                if (level == 1)
                    tagName{level} = 'x';
                else
                    tagName{level} = tagName{level-1};
                end;
            end;
            stubArr = extract_brace_string(workarr);
            mrprot = parse_loop(mrprot, stubArr, tagName, level);
            workarr = workarr(length(stubArr)+1:end);
            %if (~isempty(tagStr))
                level = level - 1;
            %end
    
        % for paramarray, special processing is required
        % for now, though, just skip them since the most important ones are
        % still in the text mrprot and it would be a lot more work to process
        % them here
        elseif (strncmpi(tagType, 'ParamArray', 10))
            if (~isempty(tagStr))
                level = level + 1;
                tagName{level} = safe_fieldname(tagStr);
            end
            stubArr = extract_brace_string(workarr);
    
    %         %skip the <default> tag, which always seems to be there
            [tagStr2, tagType, stubArr] = find_next_tag(stubArr);
            if (~strcmpi(tagStr2,'Default'))
                %fprintf('parse_xprot(): unknown ParamArray format!')
                stubArr = extract_brace_string(workarr);
            end
            % now map the fields
    %         [tagStr, tagType, stubArr] = find_next_tag(stubArr);
    
            mrprot = parse_loop(mrprot, stubArr, tagName, level);
    
            temparr = workarr;
            % walk down the string along the structure
            structpos = ['mrprot'];
            for i = 1:level
                structpos = [structpos '.' tagName{i}];
            end;
            [temparr,~,count] = walkdown(mrprot,structpos,temparr,0,0);
    
            if (length(strfind(temparr,'{')) >= count )
                % walk down the structure again, along the values of the stucture, and fill the values in
                [temparr,mrprot] = walkdown(mrprot,structpos,temparr,1,0);
            end;
    
    
            workarr = workarr(length(stubArr)+1:end);
    
            if (~isempty(tagStr))
                level = level - 1;
            end
    
        elseif (strncmpi(tagType, 'PipeService', 11))
            if (~isempty(tagStr))
                level = level + 1;
                tagName{level} = safe_fieldname(tagStr);
            end
            stubArr = extract_brace_string(workarr);
            
            %skip the <Class> tag, which always seems to be there
            [tagStr, tagType, stubArr] = find_next_tag(stubArr);
            if (~strcmpi(tagStr,'Class'))
                fprintf('parse_xprot(): unknown PipeService format!')
                stubArr = extract_brace_string(workarr);
            end
            
            mrprot = parse_loop(mrprot, stubArr, tagName, level);
            workarr = workarr(length(stubArr)+1:end);
            if (~isempty(tagStr))
                level = level - 1;
            end
            
        elseif (strncmpi(tagType, 'Pipe', 4))
            if (~isempty(tagStr))
                level = level + 1;
                tagName{level} = safe_fieldname(tagStr);
            end
            stubArr = extract_brace_string(workarr);
            mrprot = parse_loop(mrprot, stubArr, tagName, level);
            workarr = workarr(length(stubArr)+1:end);
            if (~isempty(tagStr))
                level = level - 1;
            end
            
        elseif (strncmpi(tagType, 'ParamFunctor', 12))
            if (~isempty(tagStr))
                level = level + 1;
                tagName{level} = safe_fieldname(tagStr);
            end
            stubArr = extract_brace_string(workarr);
            
            %skip the <Class> tag, which always seems to be there
            [tagStr, tagType, stubArr] = find_next_tag(stubArr);
            if (~strcmpi(tagStr,'Class'))
                fprintf('parse_xprot(): unknown PipeService format!')
                stubArr = extract_brace_string(workarr);
            end
            
            mrprot = parse_loop(mrprot, stubArr, tagName, level);
            workarr = workarr(length(stubArr)+1:end);
            if (~isempty(tagStr))
                level = level - 1;
            end
            
        % these are values that we know how to process, so do so
        elseif ((strncmpi(tagType, 'ParamBool', 9))     || ...
                (strncmpi(tagType, 'ParamLong', 9))     || ...
                (strncmpi(tagType, 'ParamDouble', 11))  || ...
                (strncmpi(tagType, 'ParamString', 11))  || ...
                (strncmpi(tagType, 'ParamChoice', 11)))
            if (isempty(tagStr))
                tagStr = 'x';
            end;
            tagName{level+1} = safe_fieldname(tagStr);
            ind = strfind(workarr,'}');
            ind2 = strfind(strtrim(workarr((ind(1)+1):end)),'{');
            if (length(ind2) == 0)
                ind2(1) = 0;
            end;
            if (ind2(1) == 1)
                value = [];
                ind = strfind(workarr,'{');
                l = 1;
                for i = 1:length(ind)
                    stubStr = extract_brace_string(workarr(ind(i):end));
                    %fprintf('%s = %s\n',tagName{level+1},stubStr);
                    temp = get_xprot_value(stubStr, tagType);  
                    if (~isempty(temp) && length(temp) == 1)
                        value(l) = temp;
                        l = l+1;
                    end;
                end;
            else
                stubStr = extract_brace_string(workarr);
                %fprintf('%s = %s\n',tagName{level+1},stubStr);
                value = get_xprot_value(stubStr, tagType);
            end;
            fields = {tagName{1:level+1}, value};
            mrprot = setfield(mrprot, fields{:});
            workarr = workarr(length(stubStr)+1:end);
    
        % we don't care about the things below, but acknowledge that we know
        % about them
        elseif (...%(strncmpi(tagType, 'ParamChoice', 11))      || ...
                (strncmpi(tagType, 'Class', 5))             || ...
                (strncmpi(tagType, 'Connection', 10))       || ...
                (strncmpi(tagType, 'Event', 5))             || ...
                (strncmpi(tagType, 'ParamFunctor', 12))     || ...
                (strncmpi(tagType, 'Method', 6))            || ...
                (strncmpi(tagType, 'ProtocolComposer', 16)) || ...
                (strncmpi(tagType, 'Dependency', 10))       || ...
                (strncmpi(tagType, 'ParamCardLayout', 15))  || ...
                (strncmpi(tagType, 'EVACardLayout', 13)))
            stubStr = extract_brace_string(workarr);
            workarr = workarr(length(stubStr)+1:end);
    
            % ignore these also
        elseif ((strncmpi(tagStr, 'Name', 4))             || ...
                (strncmpi(tagStr, 'ID', 2))               || ...
                (strncmpi(tagStr, 'Comment', 7))          || ...
                (strncmpi(tagStr, 'Label', 5))            || ...
                (strncmpi(tagStr, 'Visible', 7))          || ...
                (strncmpi(tagStr, 'Userversion', 11)))
            % skip to end of line for these
            lend = strfind(workarr, char(10));
            workarr = workarr(lend+1:end);
        elseif (strncmpi(tagStr, 'EVAStringTable', 14))
            % skip brace string for this one
            stubStr = extract_brace_string(workarr);
            workarr = workarr(length(stubStr)+1:end);
            
        % something unknown has happened if we reach this point, so throw a warning
        else
            if (~isempty(tagType))
    %             fprintf('parse_xprot(): WARNING: found unknown tag %s (%s)\n',tagStr,tagType);
            end;
        end
    
        [tagStr, tagType, workarr] = find_next_tag(workarr);
    end
end

function stvar = getQuotString(text)
    % from parse_xprot.m from E. Auerbach, CMRR, 2013
    % extracts string between double quotes, e.g. "string"
    %  also works with double-double quotes, e.g. ""string""
    
    idx = strfind(text,'"');
    
    if ( (length(idx) == 4) && (idx(1)+1 == idx(2)) && (idx(3)+1 == idx(4)) ) % double-double quotes
        stvar = text(idx(2)+1:idx(3)-1);
    elseif (length(idx) >= 2) % double quotes, or ??? just extract between first and last quotes
        stvar = text(idx(1)+1:idx(end)-1);
    else % malformed?
        stvar = text;
    end
end

function value = get_xprot_value(stubStr, tagType)
    % from parse_xprot.m from E. Auerbach, CMRR, 2013
    % stubStr contains the values of interest, but also might include
    % modifiers. we will just ignore the modifiers, which seem to always be
    % terminated by CR.
    
    [tagStr, newtagType, remStr] = find_next_tag(stubStr);
    while (~isempty(remStr)) % found a tag
        if (~isempty(newtagType)) % not a modifier tag???
            error('parse_xprot()::get_xprot_value(): ERROR: found unknown tag!')
        else
            % these are the modifier tags we know about
            if ((strncmpi(tagStr, 'Precision', 9))      || ...
                (strncmpi(tagStr, 'LimitRange', 10))    || ...
                (strncmpi(tagStr, 'MinSize', 7))        || ...
                (strncmpi(tagStr, 'MaxSize', 7))        || ...
                (strncmpi(tagStr, 'Limit', 5))          || ...
                (strncmpi(tagStr, 'Default', 7))        || ...
                (strncmpi(tagStr, 'InFile', 6))         || ...
                (strncmpi(tagStr, 'Context', 7))        || ...
                (strncmpi(tagStr, 'Dll', 3))            || ...
                (strncmpi(tagStr, 'Class', 5))          || ...
                (strncmpi(tagStr, 'Comment', 7))        || ...
                (strncmpi(tagStr, 'Label', 5))        || ...
                (strncmpi(tagStr, 'Tooltip', 7))        || ...
                (strncmpi(tagStr, 'Visible', 7))        || ...
                (strncmpi(tagStr, 'Unit', 4)))
                % acknowledge that we know about these tags
            else
    %             fprintf('parse_xprot(): WARNING: found unknown modifier %s\n',tagStr);
            end
            
            % remove the line containing the tag
            lstart = strfind(stubStr, ['<' tagStr '>']);
            %lend = strfind(stubStr(lstart+1:end), char(10));
            lend = lstart + length(tagStr) + 1;
            if (lstart > 1)
                stubStr = [stubStr(1:(lstart-1)) stubStr((lend+1):end)];
            else
                stubStr = stubStr((lend+1):end);
            end
        end
    
        [tagStr, newtagType, remStr] = find_next_tag(stubStr);
    end
    
    if (strncmpi(tagType, 'ParamChoice', 11))
        value = '';%getQuotString(stubStr);
        ok = true;
    end;
    if (strncmpi(tagType, 'ParamString', 11))
        value = getQuotString(stubStr);
        ok = true;
    else
        stubStr = strrep(stubStr, char(10), ' '); % remove newlines
        if (strncmpi(tagType, 'ParamBool', 9))
            stubStr = strrep(stubStr, '"true"', '1');
            stubStr = strrep(stubStr, '"false"', '0');
            [value,ok] = str2num(stubStr); %#ok<ST2NM>
            if (length(stubStr) > 1)
                value = false;
            else
                value = (value ~= 0);
            end;
        elseif (strncmpi(tagType, 'ParamLong', 9))
            [value,ok] = str2num(stubStr); %#ok<ST2NM>
        elseif (strncmpi(tagType, 'ParamDouble', 11))
            [value,ok] = str2num(stubStr); %#ok<ST2NM>
            ind = strfind(stubStr,'.');
            if (length(ind) > 0)
                if (length(ind) == length(value) -1) % there is an indication of bitnumber
                    value = value(2:end);
                end;
            end;
        end
    end
    
    if (~ok)
        %fprintf('WARNING: get_xprot_value failed: %s\n', stubStr);
        %disp(value);
    end
end

function [tempArr,mrprot,count] = walkdown(mrprot,structpos,tempArr,fillin,count)
    % from parse_xprot.m from E. Auerbach, CMRR, 2013
    fnames = eval(['fieldnames(' structpos ');']);   
    
    for i = 1:length(fnames)
        if (eval(['isstruct(' structpos '.' fnames{i} ')']))
            ind = strfind(tempArr,'{');
            tempArr = tempArr((ind(1)+1):end);
            [tempArr,mrprot,count] = walkdown(mrprot,[structpos '.' fnames{i}],tempArr,fillin,count);
            ind2 = strfind(tempArr,'}');
            if (isempty(ind2))
                tempArr = '';
            else
                tempArr = tempArr((ind2(1)+1):end);
            end;
            count = count +1;
        else
            ind = strfind(tempArr,'{');
            ind2 = strfind(tempArr,'}');
            if (fillin == 1)
                try
                    eval([structpos '.' fnames{i} ' = ' ( tempArr((ind(1)+1):(ind2(1)-1)) ) ';']);
                catch
                    eval([structpos '.' fnames{i} ' = 0;']);
                end;
            end;
            tempArr = tempArr((ind2(1)+1):end);
            count = count +1;
        end;
    end;
end
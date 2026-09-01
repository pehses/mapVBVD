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

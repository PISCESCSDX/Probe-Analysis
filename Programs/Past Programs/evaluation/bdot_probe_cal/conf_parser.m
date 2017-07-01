function S= conf_parser(conffile)
%function S= conf_parser(conffile)
%
% read conffile and create every found variable
% in Structure 'S' in workspace; format:
% '  varname = value   %comment '
% coments and spaces are ignored, line break with '...'
% are allowed
% if no output argument is given, the variables will
% be created in workspace of caller
%
% input:    conffile    filename
% output:   S           structure with field from file
%

    if nargin < 1
        help('conf_parser');
        return
    end

    fid= fopen(conffile, 'r');
    while 1
        wline= [];
        while 1
           tline= fgetl(fid);
           if ~ischar(tline),   break,   end
           wline= [wline, tline];
           wline= stripe_comment(wline);
           wline= stripe_spaces(wline);
           if length(wline) <= 3; break, end
           if strcmp(wline(end-2:end), '...')
               wline= wline(1:end-3);
           else
               break
           end
        end
        if ~ischar(tline),   break,   end
        
        eq_pos= findstr(wline, '=');
        if ~isempty(eq_pos)
            varname= wline(1:(eq_pos(1)-1));
            value  = wline((eq_pos(1) + 1): end);
            evalline= ['S.', varname, '= ', value, ';'];
            eval(evalline);
        end
    end
    fclose(fid);
    
    %---------------check if something was found
    anz= whos('S');
    if size(anz, 1)== 0
        warning(['No entry in file ''', conffile, ''' found.']);
    end
    
    
    
    %-----------------transfer to callers workspace
    if nargout == 0
        evalin('caller', 'global uebergabe_conf_parser;')
        global uebergabe_conf_parser;
    
        varnames= fieldnames(S);
        for i1= 1:length(varnames)
            evalline= ['uebergabe_conf_parser= S.', varnames{i1}, ';'];
            eval(evalline);
            evalline= [varnames{i1}, '= uebergabe_conf_parser;'];
            evalin('caller', evalline);
            uebergabe_conf_parser= [];
        end
        
    end
    clear global uebergabe_conf_parser
    
end

function [wline_out]= stripe_spaces(wline)
    % find positions of strings, because don't delete spaces in strings
    spos= findstr(wline, '''');
    spos= [1, spos, length(wline) + 1];
    
    wline_out= [];
    instr= false;
    i1= 0;
    while 1
        i1= i1 + 1;
        if (i1+1) > length(spos), break, end
        st= spos(i1);
        en= spos(i1 + 1) - 1;
        if ~instr
            wline_out= [wline_out, strrep(wline(st:en), ' ', '')];
        else
            wline_out= [wline_out, wline(st:en)];
        end
        instr= ~instr;
    end
end
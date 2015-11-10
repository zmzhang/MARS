currentDepth = 1; % get the supper path of the current path
currPath = fileparts(cd);% get current path
fsep = filesep;
pos_v = strfind(cd,fsep);
p = currPath(1:pos_v(length(pos_v)-currentDepth+1)-1); % -1: delete the last character '/' or '\' 
codepath=[p '\' 'code'];
path(path,codepath);
mars_main('MSRT_input','p',30,5,0.9);
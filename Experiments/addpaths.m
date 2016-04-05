mfilepath=fileparts(which(mfilename));
display(mfilepath);

addpath(genpath(fullfile(mfilepath,'\FETI')));
addpath(genpath(fullfile(mfilepath,'\FunctionsOther')));
addpath(genpath(fullfile(mfilepath,'\FunctionsOwn')));

clear mfilepath
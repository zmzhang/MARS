function mars_main(MSRT_names,project_names,mpw,w,threshold)

%      where
%
% INPUT VALUES:
%
%           mpw: max peak width (default: 30 scan points)
%             w: the window of MSWFA (default: 5 scan points)
%     threshold: eigenvalue (0.9 is the default)
% project_names: string of project names£¨defined by user)     
%    MSRT_names: string of the mat file name storting MSRT 
%

if nargin<5; threshold=0.9;end;
if nargin<4; w=5;end;
if nargin<3; mpw=50;end
if nargin<2; project_names=[];end;


% %%%%%+++++.................................reading data..............................++++%%%%% %
[filename, pathnameb, FilterIndex] = uigetfile({'*.cdf','cdf-files(*.cdf)';'*.*',...
    'All Files(*.*)'},'Pick files','MultiSelect','on'); 
full_name=fullfile(pathnameb,filename);

sam_names=cell(length(full_name),1);
for i=1:length(full_name)
    [~,name,ext] = fileparts(full_name{i});
    sam_names(i)={name};
end
eval(['save ' [project_names 'sam_names']  ' sam_names']); % save all files' names

for i=1:length(full_name)
    [~,name,ext] = fileparts(full_name{i});
    [X0,T,W]=cdftomat(full_name{i},1);   % tranfer cdf to mat
    eval(['save ' [project_names name] ' X0 T W']);    % svae data
    disp([num2str(i),'/',num2str(length(full_name)),'th of reading data']);
end

% % %%%%%+++++.....................................DECONVOLUTION.........................++++%%%%% %
load(MSRT_names)
load([project_names 'sam_names'])
sams=size(sam_names,1);
decon_names=cell(sams,1);
for i=1:sams
    name=sam_names(i);
    load(name{1})
    result = mars(X0,S,R,w,threshold,mpw);
    eval(['save ' [project_names 'd' name{1}] ' X0 T W S R result']);
    decon_names(i,1)={[project_names 'd' name{1}]};
    disp([num2str(i),'/',num2str(length(sams)),'th of resolution']);
end
eval(['save ' [project_names 'dnames'] ' decon_names']);  % save resolution results



% %%%%%+++++.....................................geting QUAL/QUAN tables.........................++++%%%%% %
load([project_names 'dnames']);
load('sam_names')
% QUAN table
A=[];
H=[];
for i=1:length(decon_names)
    name=decon_names(i);
    load(name{1});
    A=[A;result.AREAS];
    H=[H;result.HEIGS];
end
table.sam_names=sam_names;
table.A=A;table.H=H;
eval(['save ' [project_names 'QUAN_table'] ' table']);
% QUAL mass spectrums
masstomsp(project_names,S,R,W);
clear

end
%  %%%%%%%%++++++++....................................over.........................................+++++%%%%%



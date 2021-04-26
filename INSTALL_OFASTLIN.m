%--------------------------------------------------------------------------
% INSTALL_FLECCS.m
% Installs DTQP and adds the files to current MATLAB path
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Contributor: Athul K. Sundarrajan (AthulKrishnaSundarrajan on GitHub)
% Contributor: Yong Hoon Lee (yonghoonlee on GitHub)
% Contributor: Daniel R. Herber (danielrherber on GitHub)
%--------------------------------------------------------------------------
function INSTALL_OFASTLIN(varargin)

% initialize
silentflag = 0; % dont be silent

if ~isempty(varargin)
    if any(strcmpi(varargin,'silent'))
        silentflag = 1;
    end
end

% display banner text
RunSilent('DisplayBanner',silentflag)

% add contents to path
RunSilent('AddSubmissionContents(mfilename)',silentflag)

% download required web zips
RunSilent('RequiredWebZips',silentflag)

% add contents to path (files have been downloaded)
RunSilent('AddSubmissionContents(mfilename)',silentflag)

% install DTQP
RunSilent('RunThisFile',silentflag);


% open an example
if ~silentflag, OpenThisFile('Outerloop.m'); end

% close this file
RunSilent('CloseThisFile(mfilename)',silentflag)


end

function RunSilent(str,silentflag)

% if silent, capture the output
if silentflag
    O = evalc(str); %#ok<NASGU>
else
    eval(str);
end


end
function DisplayBanner

disp('--------------------------------------------------------------------------------------------')
disp('        <strong>Installing and adding required files to MATLAB path for FASTLIN</strong>     ')
disp('--------------------------------------------------------------------------------------------')
end

function AddSubmissionContents(name) %#ok<DEFNU>

disp('-> Adding submission contents to path')
disp(' ')

% turn off potential warning
warning('off','MATLAB:dispatcher:nameConflict')

% current file
fullfuncdir = which(name);

% current folder
submissiondir = fullfile(fileparts(fullfuncdir));

% add folders and subfolders to path
addpath(genpath(submissiondir))

% turn warning back on
warning('on','MATLAB:dispatcher:nameConflict')
end

function RequiredWebZips %#ok<DEFNU>

disp('-> Obtaining required web zips')

% initialize index
ind = 0;

% initialize structure
zips = struct('url','','folder','','test','');

% zip 1
ind = ind + 1; % increment
zips(ind).url = 'https://github.com/danielrherber/dt-qp-project/archive/master.zip';
zips(ind).folder = 'DTQP';
zips(ind).test = 'INSTALL_DTQP';

% zip 2
ind = ind+1; 
zips(ind).url = 'https://www.mathworks.com/matlabcentral/mlc-downloads/downloads/7c638d4c-0763-497d-8cad-dd04d42a1352/9bf146f3-da99-48f5-bae5-2cd220c3f5cc/packages/zip'; 
zips(ind).folder = 'DataFiles.zip';
zips(ind).test = 'pd_0.2_linear';
 

% obtain full function path
full_fun_path = which(mfilename('fullpath'));
outputdir = fullfile(fileparts(full_fun_path),'include');

% download and unzip
DownloadWebZips(zips,outputdir)

disp(' ')
end

function DownloadWebZips(zips,outputdir)

% store the current directory
olddir = pwd;

% create a folder for outputdir
if ~exist(outputdir, 'dir')
    mkdir(outputdir); % create the folder
else
    addpath(genpath(outputdir)); % add folders and subfolders to path
end

% change to the output directory
cd(outputdir)

% go through each zip
for k = 1:length(zips)

    % get data
    url = zips(k).url;
    folder = zips(k).folder;
    test = zips(k).test;

    % first check if the test file is in the path
    if exist(test,'file') == 0

        try
            % download zip file
            zipname = websave(folder,url);

            % save location
            if k == 2
                outputdirname = fullfile(outputdir,'DataFiles');
            else
                outputdirname = fullfile(outputdir,folder);
            end

            % create a folder utilizing name as the foldername name
            if ~exist(outputdirname, 'dir')
                mkdir(outputdirname);
            end

            % unzip the zip
            unzip(zipname,outputdirname);

            % delete the zip file
            if k ==2
                delete(folder)
            else
                delete([folder,'.zip'])
            end

            % output to the command window
            disp(['Downloaded and unzipped ',folder])

        catch % failed to download
            % output to the command window
            disp(['Failed to download ',folder])

            % remove the html file
            delete([folder,'.html'])
        end

    else
        % output to the command window
        disp(['Already available ',folder])
    end
end

% change back to the original directory
cd(olddir)
end




function OpenThisFile(name)

disp(['-> Opening ', name])

try
    % open the file
    open(name);
catch % error
    disp(['Could not open ', name])
end

disp(' ')
end
function CloseThisFile(name) %#ok<DEFNU>

disp(['-> Closing ', name])
disp(' ')

% get editor information
h = matlab.desktop.editor.getAll;

% go through all open files in the editor
for k = 1:numel(h)
    % check if this is the file
    if ~isempty(strfind(h(k).Filename,name))
        % close this file
        h(k).close
    end
end
end

function RunThisFile

% run silent 
%silentflag = 0;

% run and install DTQP
INSTALL_DTQP('silent');

end

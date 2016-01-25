function pfname = GenerateGEpfname(ImagePars)
% GenerateGEpfname - generate file names for GE image files
%
%   USAGE:
%
%   pfname = GenerateGEpfname(ImagePars)
%
%   INPUT:
%
%   ImagePars
%       Structure array that contains information about the file name
%       patterns. Must contain pname (path name), fbase (file name), fnumber(series of file
%       number), fext (file extension, i.e ge1)
%
%   OUTPUT:
% 
%   pfname
%       cell structure with the full file names
% 

numimages   = length(ImagePars.fnumber);
for i = 1:1:numimages
    fname   = sprintf([ImagePars.fbase, '%05d.%s', ], ImagePars.fnumber(i), ImagePars.fext);
    pfname{i,1} = fullfile(ImagePars.pname, fname);
end

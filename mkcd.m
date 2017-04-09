function [success] = mkcd(directory)
%MKCD('directory') makes directory and cds to it
mkdir(directory)
cd(directory)
success =true;

end


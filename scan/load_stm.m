function [stm]=load_stm(filename)
% Load stimulator timestamps from trial.stm file.

fid = fopen([filename '.stm'],'r','ieee-be');        % Open tet file. 'ieee-be' string is machine format, which needs to be 'big endian' to read times correctly. Default format doesnt work on my PC.
if fid==-1
    stm=[];
    return
end

hdr = fread(fid,400,'int8');
hdr=char(abs(hdr))';

timebase_index = (findstr('timebase ',hdr))+8;
timebase = str2double(strtok(hdr(timebase_index:end)));

nStm_index = (findstr('num_stm_samples',hdr))+15;
nStm = str2double(strtok(hdr(nStm_index:end)));

ds = (findstr('data_start',hdr))+9;        % Look for data start marker

fseek(fid,ds,'bof');

stm=fread(fid,nStm,'uint32');

stm=stm./timebase;

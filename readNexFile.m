function nex = readNexFile(fname, verbose)
% readNexFile - Reads a NEX file (Nex Technologies), and output the data in
% a MATLAB structure. Wrapper script for readNexFileC.
%
% USAGE:
%   nex = readNexFile('filename', verbose)
%
% readNexFile reads a NEX file (Nex Technologies), and output the data in
% a MATLAB structure. This is actually a wrapper script for readNexFileC.
% This wrapper takes the output from readNexFileC, and converts the data
% into a more user friendly format. In particular, it converts all
% variables to doubles (from either int32 or int16). It also converts raw
% ticks into seconds, and raw A/D values into millivolts (mV). See
% readNexFileC.c fo
%
% INPUT:
%   filename - Name of the NEX file to read.
%   verbose  - Enable verbose output (optional, default 'false')
%
% OUTPUT:
%   nex - A structure containing .nex file data
%   nex.version - File version
%   nex.comment - File comment
%   nex.freq    - Frequency (to convert from ticks to seconds)
%   nex.tbeg    - Beginning of recording session (in seconds)
%   nex.tend    - End of recording session (in seconds)
%
%   nex.neurons - Cell array of neuron structures
%       neuron.name       - Name of neuron variable
%       neuron.timestamps - Array of neuron timestamps (in seconds)
%
%   nex.events - Cell array of event structures
%       event.name       - Name of event variable
%       event.timestamps - Array of event timestamps (in seconds)
%
%   nex.intervals - Cell array of interval structures
%       interval.name      - Name of interval variable
%       interval.intStarts - Array of interval starts (in seconds)
%       interval.intEnds   - Array of interval ends (in seconds)
%
%   nex.waves - Cell array of waveform structures
%       wave.name        - Name of waveform variable
%       wave.WFrequency  - A/D frequency for waveform data points
%       wave.ADtoMV      - Conversion from A/D values to mV
%       wave.NPointsWave - Number of data points in each wave
%       wave.MVOffset    - Offset of A/D values
%                          (mv = raw*ADtoMV+MVOffset)
%       wave.timestamps  - Array of waveform timestamps (in seconds)
%       wave.waveforms   - Matrix of waveforms (in mV),
%                          each column vector is a wave.
%
%   nex.popvectors - Cell array of population vector structures
%       popvector.name    - Name of population vector variable
%       popvector.weights - Array of population vector weights
%
%   nex.contvars - Cell array of continuous variable structures
%       contvar.name         - Name of continuous variable
%       contvar.ADFrequency  - A/D frequency for data points
%       contvar.ADtoMV       - Conversion from A/D values to mV
%       contvar.MVOffset     - Offset of A/D values
%                              (mv =raw*ADtoMV+MVOffset)
%
%       continuous (A/D) data come in fragments. Each fragment has a
%       timestamp and an index of the A/D data points in data array. The
%       timestamp corresponds to the time of recording of the first A/D
%       value in this fragment.
%
%       contvar.timestamps     - Array of timestamps (in seconds)
%       contvar.fragmentStarts - Array of start indexes
%       contvar.data           - Array of data points (in mV)
%
%   nex.markers - Cell array of marker structures
%       marker.name              - Name of marker variable
%       marker.length            - Number of characters per marker value
%       marker.timestamps        - Array of marker timestamps (in seconds)
%       marker.values            - Array of marker value structures
%           marker.value.name    - Name of marker value
%           marker.value.strings - Array of marker value strings
%
% See also writeNexFile readNexFileC writeNexFileC
%
% AUTHOR: Benjamin Kraus (bkraus@bu.edu, ben@benkraus.com)
% Copyright (c) 2008-2010, Benjamin Kraus
% $Id: readNexFile.m 1104 2010-01-16 22:36:42Z bkraus $

nex = [];

if (nargin == 0 || isempty(fname))
    FilterSpec = {'*.nex', 'NeuroExplorer File (*.nex)';
        '*', 'All Files'};
    [fname, pathname] = uigetfile(FilterSpec, 'Select a NeuroExplorer file');
    if(fname == 0); return; end
    fname = strcat(pathname, fname);
end

if(nargin <= 1); verbose = false; end

input = readNexFileC(fname, verbose);

nex = input;
freq = input.freq;
nex.version = double(input.version);
nex.comment = deblank(input.comment);
nex.tbeg = double(input.tbeg)./freq;
nex.tend = double(input.tend)./freq;

if(isfield(input,'neurons'))
    for i = 1:size(input.neurons,1)
        if(isempty(nex.neurons{i,1})); break; end
        nex.neurons{i,1}.version = double(input.neurons{i,1}.version);
        nex.neurons{i,1}.wireNumber = double(input.neurons{i,1}.wireNumber);
        nex.neurons{i,1}.unitNumber = double(input.neurons{i,1}.unitNumber);
        nex.neurons{i,1}.gain = double(input.neurons{i,1}.gain);
        nex.neurons{i,1}.filter = double(input.neurons{i,1}.filter);
        nex.neurons{i,1}.timestamps = double(input.neurons{i,1}.timestamps)./freq;
    end
    if(isempty(input.neurons)); nex = rmfield(nex,'neurons'); end
end

if(isfield(input,'events'))
    for i = 1:size(input.events,1)
        if(isempty(nex.events{i,1})); break; end
        nex.events{i,1}.version = double(input.events{i,1}.version);
        nex.events{i,1}.timestamps = double(input.events{i,1}.timestamps)./freq;
    end
    if(isempty(input.events)); nex = rmfield(nex,'events'); end
end

if(isfield(input,'intervals'))
    for i = 1:size(input.intervals,1)
        if(isempty(nex.intervals{i,1})); break; end
        nex.intervals{i,1}.version = double(input.intervals{i,1}.version);
        nex.intervals{i,1}.intStarts = double(input.intervals{i,1}.intStarts)./freq;
        nex.intervals{i,1}.intEnds = double(input.intervals{i,1}.intEnds)./freq;
    end
    if(isempty(input.intervals)); nex = rmfield(nex,'intervals'); end
end

if(isfield(input,'waves'))
    for i = 1:size(input.waves,1)
        if(isempty(nex.waves{i,1})); break; end
        ADtoMV = input.waves{i,1}.ADtoMV;
        MVOffset = input.waves{i,1}.MVOffset;
        nex.waves{i,1}.version = double(input.waves{i,1}.version);
        nex.waves{i,1}.wireNumber = double(input.waves{i,1}.wireNumber);
        nex.waves{i,1}.unitNumber = double(input.waves{i,1}.unitNumber);
        nex.waves{i,1}.NPointsWave = double(input.waves{i,1}.NPointsWave);
        nex.waves{i,1}.timestamps = double(input.waves{i,1}.timestamps)./freq;
        nex.waves{i,1}.waveforms  = double(input.waves{i,1}.waveforms).*ADtoMV + MVOffset;
    end
    if(isempty(input.waves)); nex = rmfield(nex,'waves'); end
end

if(isfield(input,'popvectors'))
    for i = 1:size(input.popvectors,1)
        if(isempty(nex.popvectors{i,1})); break; end
        nex.popvectors{i,1}.version = double(input.popvectors{i,1}.version);
    end
    if(isempty(input.popvectors)); nex = rmfield(nex,'popvectors'); end
end

if(isfield(input,'contvars'))
    for i = 1:size(input.contvars,1)
        if(isempty(nex.contvars{i,1})); break; end
        ADtoMV = input.contvars{i,1}.ADtoMV;
        MVOffset = input.contvars{i,1}.MVOffset;
        nex.contvars{i,1}.version = double(input.contvars{i,1}.version);
        nex.contvars{i,1}.timestamps = double(input.contvars{i,1}.timestamps)./freq;
        nex.contvars{i,1}.fragmentStarts = double(input.contvars{i,1}.fragmentStarts)+1;
        nex.contvars{i,1}.data = double(input.contvars{i,1}.data).*ADtoMV + MVOffset;
    end
    if(isempty(input.contvars)); nex = rmfield(nex,'contvars'); end
end

if(isfield(input,'markers'))
    for i = 1:size(input.markers,1)
        if(isempty(nex.markers{i,1})); break; end
        nex.markers{i,1}.version = double(input.markers{i,1}.version);
        nex.markers{i,1}.length = double(input.markers{i,1}.length);
        nex.markers{i,1}.timestamps = double(input.markers{i,1}.timestamps)./freq;
        
        for j= 1:size(input.markers{i,1}.values,1);
            nex.markers{i,1}.values{j,1}.name = deblank(input.markers{i,1}.values{j,1}.name);
            nex.markers{i,1}.values{j,1}.strings = deblank(cellstr(input.markers{i,1}.values{j,1}.strings));
        end
    end
    if(isempty(input.markers)); nex = rmfield(nex,'markers'); end
end
end

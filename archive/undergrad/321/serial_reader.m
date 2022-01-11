function [t, data] = serial_reader(numDataPoints)
% Copyright 2014 The MathWorks, Inc.
% with some modificaitons by Eric Mehiel

%% Create serial object for Arduino
if (~isempty(instrfind))
    fclose(instrfind); % close all ports to start, just to be sure...
    delete(instrfind);
end

s = serial('COM11');  % change the COM Port number as needed

%% Connect the serial port to Arduino

try
    fopen(s);
catch err
    fclose(instrfind);
    delete(instrfind);
    error('Make sure you select the correct COM Port where the Arduino is connected.');
end

%% Read  the data from Arduino

ii = 0; % counter

data = zeros(numDataPoints,5);
t = zeros(numDataPoints,1);


timer = tic; % Start timer
while ii < numDataPoints
    ii = ii + 1;
    % Read buffer data
    disp('Reading Serial Data');
    data(ii,:) = fscanf(s, '%e %e %e %e %e')'; % Change format string as needed
    % Read time stamp
    t(ii) = toc(timer);
end
fclose(s);
delete(s);
clear s;
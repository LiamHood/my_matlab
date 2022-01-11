function [t, data] = serial_reader(numDataPoints)
% Copyright 2014 The MathWorks, Inc.
% with some modificaitons by Eric Mehiel

%% Create serial object for Arduino
if (~isempty(instrfind))
    fclose(instrfind); % close all ports to start, just to be sure...
    delete(instrfind);
end

s = serial('COM6');  % change the COM Port number as needed

%% Connect the serial port to Arduino

try
    fopen(s);
catch err
    fclose(instrfind);
    delete(instrfind);
    error('Make sure you select the correct COM Port where the Arduino is connected.');
end

%% Read  the data from Arduino

i = 0; % counter
%numDataPoints = 10;
data = zeros(numDataPoints,2);
t = zeros(numDataPoints,1);

tempData = fscanf(s, '%e'); % Clear the serial buffer

timer = tic; % Start timer
while i < numDataPoints
    i = i + 1;
    % Read buffer data
    disp('Reading Serial Data');
    data(i,:) = fscanf(s, '%e %e')'; % Change format string as needed
    % Read time stamp
    t(i) = toc(timer);
end
fclose(s);
delete(s);
clear s;
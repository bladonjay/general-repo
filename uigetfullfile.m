function [fullfile] = uigetfullfile(varargin)
% function [fullfile] = uigetfullfile();
% gets the full file out of uigetfile


[a,b]=uigetfile();
fullfile=[b a];


end


% commaDelimitedToCell.m
%
%        $Id: commaDelimitedToCell.m 666 2007-10-05 14:15:21Z justin $
%      usage: cellarray = commaDelimitedToCell(str)
%         by: justin gardner
%       date: 10/05/07
%    purpose: convert a comma delimited string to a cell array
%
function out = commaDelimitedToCell(in)

% check arguments
if ~any(nargin == [1])
  help commaDelimitedToCell
  return
end

out = {};
while ~isempty(in)
  [out{end+1} in] = strtok(in,',');
end

% cdanat.m
%
%        $Id: cdanat.m 1166 2008-05-28 15:26:29Z justin $
%      usage: cdanat()
%         by: justin gardner
%       date: 06/19/07
%    purpose: switch to anat directory
%
function retval = cdanat()

% check arguments
if ~any(nargin == [0])
  help cdanat.m
  return
end

volDir = mrGetPref('volumeDirectory');
if isdir(volDir)
  cd(volDir);
else
  disp(sprintf('(cdanat) Could not find directory %s',volDir));
end


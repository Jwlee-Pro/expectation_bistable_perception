% mrParamsSet.m
%
%        $Id: mrParamsSet.m 1289 2008-08-25 20:15:15Z justin $
%      usage: mrParamsSet(params)
%         by: justin gardner
%       date: 03/18/08
%    purpose: pass in params and this will reset the fields in an
%             open mrParamsDialog, useful for when you have a
%             non-modal dialog. Note that you do not have to havee
%             all the fields in the params structure, only the ones
%             you want to set
% 
%
function retval = mrParamsSet(params)

% check arguments
if ~any(nargin == [1])
  help mrParamsSet
  return
end

global gParams;


% go through each one of the passed in parameters
paramFields = fieldnames(params);
for fnum = 1:length(paramFields)
  % look for the parameter in varinfo
  match = 0;
  for vnum = 1:length(gParams.varinfo)
    if strcmp(paramFields{fnum},gParams.varinfo{vnum}.name)
      match = vnum;
    end
  end
  % found the match, go set the field
  if match
    % numeric
    if strcmp(gParams.varinfo{match}.type,'numeric')
      % check if this is a group
      if isfield(gParams.varinfo{match},'group')
	groupVar = gParams.varinfo{match}.group;
	% keep values for each one of the group
	for i = 1:length(params.(paramFields{fnum}))
	  gParams.varinfo{match}.allValues{i} = params.(paramFields{fnum})(i);
	end
	% and set the current value in the gui
	currentVal = params.(paramFields{fnum})(params.(groupVar));
	set(gParams.ui.varentry{match},'String',num2str(currentVal));
      
      else
	% simple numeric
	currentVal = params.(paramFields{fnum});
	set(gParams.ui.varentry{match},'String',num2str(currentVal));
      end
      % check for minmax, so that we can turn arrows off and on appropriately
      if isfield(gParams.varinfo{match},'minmax')
	minmax = gParams.varinfo{match}.minmax;
	if isfield(gParams.ui,'incdec')
	  incdecUI = gParams.ui.incdec{match};
	  if ~isempty(incdecUI)
	    if currentVal == min(minmax)
	      set(incdecUI(1),'Enable','off');
	    else
	      set(incdecUI(1),'Enable','on');
	    end
	    if currentVal == max(minmax)
	      set(incdecUI(2),'Enable','off');
	    else
	      set(incdecUI(2),'Enable','on');
	    end
	  end
	end
      end	  
    % checkbox
    elseif strcmp(gParams.varinfo{match}.type,'checkbox')
      set(gParams.ui.varentry{match},'Value',params.(paramFields{fnum}));
      % array
    elseif strcmp(gParams.varinfo{match}.type,'array')
      if isequal(size(gParams.varinfo{match}.value),size(params.(paramFields{fnum})))
	for matrixX = 1:size(params.(paramFields{fnum}),1)
	  for matrixY = 1:size(params.(paramFields{fnum}),2)
	    set(gParams.ui.varentry{match}(matrixX,matrixY),'String',num2str(params.(paramFields{fnum})(matrixX,matrixY)));
	  end
	end
      else
	disp(sprintf('(mrParamsSet) Array size of variable %s does not match',paramsFields{fnum}));
      end
    % popupmenu
    elseif strcmp(gParams.varinfo{match}.type,'popupmenu')
      % check if this is a group
      if isfield(gParams.varinfo{match},'group')
	% Don't do anything for a group parameter--fix later if you need this
      else
	value = find(strcmp(params.(paramFields{fnum}),gParams.varinfo{match}.value));
	if ~isempty(value)
	  set(gParams.ui.varentry{match},'Value',value);
	else
	  disp(sprintf('(mrParamsSet) %s is not a valid option for variable %s',params.(paramFields{fnum}),paramFields{fnum}));
	end
      end
    % push button. Don't do anything
    elseif strcmp(gParams.varinfo{match}.type,'pushbutton')
    % string
    elseif strcmp(gParams.varinfo{match}.type,'string')
      % check if we are group'd
      if isfield(gParams.varinfo{match},'group')
	groupVar = gParams.varinfo{match}.group;
	% keep values for each one of the group
	for i = 1:length(params.(paramFields{fnum}))
	  gParams.varinfo{match}.allValues{i} = params.(paramFields{fnum}){i};
	end
	stringValue = params.(paramFields{fnum}){params.(groupVar)};
      % check if we are contingentOn
      elseif isfield(gParams.varinfo{match},'contingentOn')
	contingentValue = params.(gParams.varinfo{gParams.varinfo{match}.contingentOn}.name);
	stringValue = gParams.varinfo{match}.allValues{contingentValue};
      else
	stringValue = params.(paramFields{fnum});
      end
      set(gParams.ui.varentry{match},'String',stringValue);
    % unimplemented type
    else
      disp(sprintf('(mrParamsSet) Setting of type %s not implemented yet',gParams.varinfo{match}.type));
    end
  else
    if ~strcmp(paramFields{fnum},'paramInfo')
      disp(sprintf('(mrParamsSet) Could not find var %s',paramFields{fnum}));
    end
  end
end



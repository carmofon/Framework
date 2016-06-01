function varargout = problem_path(varargin)
global mypath;
% PROBLEM_PATH Application M-file for problem_path.fig
%    FIG = PROBLEM_PATH launch problem_path GUI.
%    PROBLEM_PATH('callback_name', ...) invoke the named callback.

% Last Modified by GUIDE v2.0 30-Mar-2003 22:09:34

if nargin == 0  % LAUNCH GUI

	fig = openfig(mfilename,'reuse');

	% Generate a structure of handles to pass to callbacks, and store it. 
	handles = guihandles(fig);
	guidata(fig, handles);

    % Load the listbox
    mypath='C:';
    update_listbox(handles);
    update_popup(handles);
    
	if nargout > 0
		varargout{1} = fig;
	end

elseif ischar(varargin{1}) % INVOKE NAMED SUBFUNCTION OR CALLBACK

	try
		[varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard
	catch
		disp(lasterr);
	end

end


%| ABOUT CALLBACKS:
%| GUIDE automatically appends subfunction prototypes to this file, and 
%| sets objects' callback properties to call them through the FEVAL 
%| switchyard above. This comment describes that mechanism.
%|
%| Each callback subfunction declaration has the following form:
%| <SUBFUNCTION_NAME>(H, EVENTDATA, HANDLES, VARARGIN)
%|
%| The subfunction name is composed using the object's Tag and the 
%| callback type separated by '_', e.g. 'slider2_Callback',
%| 'figure1_CloseRequestFcn', 'axis1_ButtondownFcn'.
%|
%| H is the callback object's handle (obtained using GCBO).
%|
%| EVENTDATA is empty, but reserved for future use.
%|
%| HANDLES is a structure containing handles of components in GUI using
%| tags as fieldnames, e.g. handles.figure1, handles.slider2. This
%| structure is created at GUI startup using GUIHANDLES and stored in
%| the figure's application data using GUIDATA. A copy of the structure
%| is passed to each callback.  You can store additional information in
%| this structure at GUI startup, and you can change the structure
%| during callbacks.  Call guidata(h, handles) after changing your
%| copy to replace the stored original so that subsequent callbacks see
%| the updates. Type "help guihandles" and "help guidata" for more
%| information.
%|
%| VARARGIN contains any extra arguments you have passed to the
%| callback. Specify the extra arguments by editing the callback
%| property in the inspector. By default, GUIDE sets the property to:
%| <MFILENAME>('<SUBFUNCTION_NAME>', gcbo, [], guidata(gcbo))
%| Add any extra arguments after the last argument, before the final
%| closing parenthesis.

% --------------------------------------------------------------------
function varargout = txtInfo_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.txtInfo.

% --------------------------------------------------------------------
function varargout = Lista_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.Lista.
global mypath;
global file_path;

index_selected = get(handles.Lista,'Value');

if index_selected == 2
    end_backslash = 0;
    begin_point = 0;
    for i = 1:length(mypath)
        if strcmp(mypath(i),'\')
            end_backslash = i;
        end
    end
    mypath = mypath(1 : end_backslash - 1);
else
    begin_point = 0;
    adauga = strcat('\',file_path{index_selected});
    for j = 1:length(adauga)
        if strcmp(adauga(j),'.')
            begin_point = j;
            break
        end
    end
    if begin_point == 0
        mypath = [mypath adauga];
    end
end

set(handles.txtInfo,'String',mypath);
if begin_point == 0
    update_listbox(handles);
end

% --------------------------------------------------------------------
function varargout = btnApply_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.btnApply.
global mypath;
mypath = get(handles.txtInfo,'String');
close(handles.browser);

% --------------------------------------------------------------------
function varargout = btnCancel_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.btnCancel.
close(handles.browser);

% --------------------------------------------------------------------
function update_listbox(handles)
global mypath;
global file_path;

files = dir(mypath);
k=0;
x = []; begin_point = 0; j = 0;
for i = 1:max(size(files))
    if files(i).isdir
    k=k+1;
    file_path{k} = files(i).name;
    for j = 1:length(files(i).name)
        filename = files(i).name;
        if strcmp(filename(j),'.')
            begin_point = j;
        end
    end
    if (i == 1) | (i == 2)
        x{k} = ['              ' files(i).name '\'];
    elseif (i ~= 1) & (i ~= 2) & (begin_point == 0)
        x{k} = ['<DIR>     ',files(i).name];
    else
        x{k} = ['              ',files(i).name];
    end
    begin_point = 0;
end
    
end

set(handles.Lista,'String',x,'Value',1);
set(handles.txtInfo,'String',mypath);




% --------------------------------------------------------------------
function varargout = popupmenu1_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.popupmenu1.
global mypath;

val=get(handles.popupmenu1,'Value');
select_list=get(handles.popupmenu1,'String');
mypath=select_list{val}
set(handles.txtInfo,'String',mypath);
update_listbox(handles);

% --------------------------------------------------------------------
function update_popup(handles)
global mypath;

ret = {};
   
startletter = 'c';
for i = startletter:'z'
    if exist([i ':\']) == 7
        ret{end+1} = [i ':'];
    end
end
set(handles.popupmenu1,'String',ret);
mypath=ret{1};


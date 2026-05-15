function [constriction,dilation, done]=select_capillaries(c,p,valc,Name,Folder,constriction,dilation)

global r ht hl name done f Folder2
done=0;
name=Name;
Folder2=Folder;
% unrestricted radius
r=zeros(1,size(c,1));
% set existing values and treat passed constriction/dilation when caller provided them
if nargin >= 6
    if ~isempty(constriction)
    r(constriction) = -1;
    end
end
if nargin >= 7
    if ~isempty(dilation)
    r(dilation) = 1;
    end
end

% open new figure
f=figure;
clf
hold on
uicontrol(f,'Style','pushbutton','String','Done',...
    'Position',[20 20 80 40],'Callback',@done_data);
uicontrol(f,'Style','pushbutton','String','Load',...
    'Position',[20 70 80 40],'Callback',@load_data);
ht=uicontrol(f,'Style','text','Position',[20 180 80 40],'horizontalalignment','left');
uicontrol(f,'Style','text','Position',[20 240 120 80],'horizontalalignment','left',...
    'string',{'Left click on vessel: constrict','','Right click on vessel: dilate'});

% plot data
x=[p(c(:,1),1) p(c(:,2),1)];
y=[p(c(:,1),2) p(c(:,2),2)];
hl=zeros(1,length(r));
valc=reshape(valc,1,[]);
for i=valc
    % plot variable vessel and store index
    hl(i)=plot(x(i,:),y(i,:));
    set(hl(i),'ButtonDownFcn',@mbf,'Userdata',i)
end
for i=setdiff(1:size(c,1),valc)
    % plot fixed vessel
    plot(x(i,:),y(i,:),'color',[0.4 0.4 0.4],'linewidth',4);
end
axis off equal
set(gca,'Ydir','reverse')
xl=[min(x(:)) max(x(:))];
yl=[min(y(:)) max(y(:))];
xlim(xl+0.02*[-1 1]*diff(xl));
ylim(yl+0.02*[-1 1]*diff(yl));
updateplot

% wait till done button is pressed or window is closed
uiwait(gcf)
constriction=find(r==-1);
dilation=find(r==1);
if ishandle(f)
    close(f)
end



function updateplot
global r ht hl
col=[0 0 0.5;0 0 1;0.5 0.5 1];
set(ht,'string',sprintf('# constricted: %g\n\n# dilated: %g',sum(r==-1),sum(r==1)))
for i=1:length(r)
    if hl(i)
        set(hl(i),'linewidth',r(i)*3+4,'color',col(r(i)+2,:))
    end
end
shg


function mbf(hObject,eventdata)
global r

i=get(hObject,'Userdata');
% r0=r;
switch get(ancestor(hObject,'figure'),'SelectionType')
    case 'normal' % left click: constriction
        r(i)=max(r(i)-1,-1);
    case 'alt' % right click: dilation
        r(i)=min(r(i)+1,1);
end
% check maximum number of changes
% if sum(r~=0)>50
%     r=r0;
% end
updateplot


function load_data(hObject,eventdata)
global r
[FileName,PathName] = uigetfile('*.mat','Select the MAT file to load');
load([PathName,FileName])
r(:)=0;
r(constriction)=-1;
r(dilation)=1;
updateplot


function done_data(hObject,eventdata)
global r name done Folder2
Path=cd;
done=1;
constriction=find(r==-1);
dilation=find(r==1);
save([Folder2 char(name) '_Setting.mat'],'constriction','dilation')
uiresume(gcbf)



function fig_h = zgrid_hires( Handle, ReqStr )

% Hopkins function:  zgrid_hires
%
% Draws a high-resolution z-plane damping/natural-frequency grid,
% into a specified axes or figure. 
% If no axes is specified, draws a hi-res grid in the current axes.
% If no axes exists in the current figure, creates an axes first.
% If there is no current figure, creates a new figure and axes, 
% then draws a high-res grid into the new axes.
%
% Format:  figure_number = zgrid_hires( < figure_number > )
%
% in which the optional input argument
%     figure_number   <--  figure number (or axes handle)
%
% (c)2004-2014, Mark A. Hopkins
% ===================================================================

global PZG

fig_h = [];
if ~nargin
  if ~isempty(gcbf)
    ax_h = findobj( gcbf,'type','axes');
    for k = numel(ax_h):-1:1
      if strcmp( get(ax_h(k),'tag'),'legend')
        ax_h(k) = [];
      end
    end
    Handle = ax_h;
  else
    Handle = [];
  end
  if isempty(Handle)
    Handle = gcf;
  end
end

if numel(Handle) ~= 1
  disp('ZGRID_HIRES Error:  First input arg must be scalar figure handle.')
  return
end  
if ~ishandle(Handle) ...
  ||( ~strcmp( get(Handle,'type'),'axes') ...
     && ~strcmp( get(Handle,'type'),'figure') )
  disp('ZGRID_HIRES Error:  Input arg is not a figure or axes handle.')
  return
elseif strcmp( get(Handle,'type'),'figure')
  ax_h = findobj( Handle,'type','axes');
  if isempty(ax_h)
    Handle = axes('parent', Handle );
  else
    for k = numel(ax_h):-1:1
      if strcmp( get(ax_h(k),'tag'),'legend')
        ax_h(k) = [];
      end
    end
    Handle = ax_h;
  end
  
  if length(Handle) > 1
    disp('ZGRID_HIRES Error:  Specified figure has more than one axes.')
    return
  elseif isempty(Handle)
    if nargin
      disp('ZGRID_HIRES Error:  Specified figure has no axes.')
    else
      disp('ZGRID_HIRES Error:  Current figure has no axes.')
    end
    return
  end  
end

if strcmp( get(Handle,'type'),'figure');
  fig_h = Handle;
elseif strcmp( get(Handle,'type'),'axes');
  fig_h = get( Handle,'parent');
end

if nargin == 2
  if ischar(ReqStr)
    return
  end
  if ~ischar(ReqStr)
    if isnumeric(ReqStr)
      Fs = ReqStr(1);
      ReqStr = 'activate mouse';
      set( Handle,'userdata', Fs )
    else
      return
    end
  else
    Fs = get( Handle,'userdata');
    if ~isnumeric(Fs)
      Fs = [];
    end
    if isempty(Fs)
      if strcmp( get(Handle,'type'),'figure')
        Fs = get( get( Handle,'currentaxes'),'userdata');
      elseif strcmp( get(Handle,'type'),'axes')
        Fs = get( get(Handle,'parent'),'userdata');
      end
      if ~isnumeric(Fs)
        Fs = [];
      elseif numel(Fs) > 1
        Fs = [];
      end
    end
    if isempty(Fs)
      if length(PZG) == 2
        if isfield( PZG(2),'Ts')
          Fs = 1/PZG(2).Ts;
        end
      end
      if isempty(Fs)
        return
      end
    end
  end
  YLim = get( Handle,'ylim');
  XLim = get( Handle,'xlim');
  switch ReqStr
  case 'mouse motion'
    temp = get( Handle,'currentpoint');
    if ( temp(1,1) > XLim(1) ) && ( temp(1,1) < XLim(2) ) ...
      &&( temp(1,2) > YLim(1) ) && ( temp(1,2) < YLim(2) )
      CurrPt = temp(1,1) + 1i*temp(1,2);
      log_loc = log(CurrPt);
      DampPct = -100*sign(real(log_loc)) ...
                .* 1./abs( 1 + 1i*imag(log_loc)./real(log_loc) );
      NatFrqHz = abs(log(CurrPt))*Fs/2/pi;
      DispStr = {['Zeta = ' num2str(DampPct) '%']; ...
                 ['Wn = ' num2str(NatFrqHz) ' Hz']};
      textH = findobj( Handle,'type','text');
      if isempty(textH)
        if strcmp( get(Handle,'type'),'figure')
          AxHandle = get( Handle,'currentaxes');
          textH = text( 0.6*XLim(1)+0.4*XLim(2), 0.8*YLim(1)+0.2*YLim(2), DispStr,'parent', AxHandle );
        elseif strcmp( get(Handle,'type'),'axes')
          textH = text( 0.6*XLim(1)+0.4*XLim(2), 0.8*YLim(1)+0.2*YLim(2), DispStr,'parent', Handle );
        end
      else
        if length(textH) > 1
          delete(textH(2:end))
          textH = textH(1);
        end
        set( textH,'string', DispStr )
        textP = get( textH,'position');
        if ( textP(1) ~= (0.6*XLim(1)+0.4*XLim(2)) ) ...
          ||( textP(2) ~= (0.8*YLim(1)+0.2*YLim(2)) )
          set( textH,'position',[ 0.6*XLim(1)+0.4*XLim(2), ...
                                  0.8*YLim(1)+0.2*YLim(2), 0 ] )
        end
      end
      set( textH,'fontweight','bold','color','r','fontsize',14)
    else
      textH = findobj( Handle,'type','text');
      if ~isempty(textH)
        if length(textH) > 1
          delete(textH(2:end))
          textH(2:end) = [];
        end
        if ~sum(get(get(textH,'parent'),'color'))
          set( textH,'color',[ 1 1 1] )
        end
        textS = get(textH,'string');
        if ~isempty(textS)
          set(textH,'string','')
        end
      end
    end
    pause(0.05)
    return
  case 'mouse button'
    return
  case 'activate mouse'
    FigH = get( Handle,'parent');
    if isempty(get(FigH,'windowbuttondownfcn')) ...
      && isempty(get(FigH,'windowbuttonmotionfcn'))
      set( FigH,'windowbuttondownfcn', ...
                   ['zoom(gcbf,''down'');', ...
                    'zgrid_hires(get(gcbf,''currentaxes''),''mouse button'');'], ...
                'windowbuttonmotionfcn', ...
                   'zgrid_hires(get(gcbf,''currentaxes''),''mouse motion'');')
    end
  end
end

%htemp = axes(Handle); %#ok<NASGU>
delH = findobj(Handle,'tag','HiRes Zgrid');
delete(delH)

Zetas = (0.1:0.1:0.9)';
Zetas = [ 0; 0.01; 0.02; 0.05; Zetas];
normPt = -Zetas + 1i*sqrt(1-Zetas.^2);
LogBase = [ (0:0.0001:0.0099) ...
            logspace(-2,0,250) ...
            logspace(0,log10(pi),150)]';
ZLog = zeros(size(LogBase));
ZLog(:,1)=exp(normPt(1)*LogBase);
for Ck = 2:length(Zetas)
  ZLog(:,Ck)=exp(normPt(Ck)*LogBase/imag(normPt(Ck)));
end

UC = ZLog(:,1);
UC = [ 1i; 0; UC; conj(UC); 0; -1i ];

ZLog(:,1) = [];
ZLog = [ ZLog; NaN*ones(1,size(ZLog,2)) ];
ZLog = [ ZLog(:); conj(ZLog(:)) ];

set( Handle,'nextplot','add')
% plot( ZLog, ...
%      'color',[0 0 0], ...
%      'parent', Handle, ...
%      'linestyle',':', ...
%      'linewidth', 1, ...
%      'tag','HiRes Zgrid' );
plot( UC, ...
     'color',[0 0.7 0.7], ...
     'linestyle','-', ...
     'linewidth', 1, ...
     'parent', Handle, ...
     'tag','HiRes Zgrid' );


% set(H(4),'linewidth', 1,'linestyle','-','color',[0 0.7 0.7]);
% H=plot(conj(ZLog),'k:','parent', Handle );
% set(H,'linewidth', 1,'tag','HiRes Zgrid');
% set(H(4),'linewidth', 1,'linestyle','-','color',[0 0.7 0.7]);
 
LogBase = logspace(-3,log10(pi/2),300)';
normPt = pi*( -cos(LogBase)+1i*sin(LogBase) );
WLog = exp( 0.1*normPt );
for Ck = 2:10
  WLog = [WLog exp( Ck/10*normPt )];     %#ok<AGROW>
end
WLog = [ WLog; NaN*ones(1,size(WLog,2)) ];
WLog = [ WLog(:); conj(WLog(:)) ];

% H = plot( WLog(:),'k:','parent', Handle );
% set( H,'linewidth', 1,'tag','HiRes Zgrid');
plot( [ ZLog; WLog ], ...
     'color',[0 0 0], ...
     'parent', Handle, ...
     'linestyle',':', ...
     'linewidth', 1, ...
     'tag','HiRes Zgrid' );

% H = plot(conj(WLog),'k:','parent', Handle );
% set( H,'linewidth', 1,'tag','HiRes Zgrid');

set( Handle,'xgrid','on','ygrid','on')

return


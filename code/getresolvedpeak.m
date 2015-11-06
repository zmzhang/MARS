function getresolvedpeak(project_names,sam,com,f)

load([project_names 'dnames']);
load(decon_names{sam});
tzs=result.TZ;
tz=tzs(1,:,com);
cs=result.C;
c=cs{com};
s=S(com,:);
if max(tz)==0
    disp('no such peak existed')
elseif max(c)==0
    disp('peaks can not deconvolute')
    minrz=min(min(tz));
    maxrz=max(max(tz));
    tzlab1=(minrz+min(T)/f-1)*f;
    tzlab2=(maxrz+min(T)/f-1)*f;
    yy=X0(minrz:maxrz,:);
    maxmax=1.1*max(max(yy));
    plot(tzlab1:f:tzlab2,yy);
    axis([tzlab1 tzlab2 0 maxmax]);
    xlabel('elution time (s)')
    ylabel('intensity')
    title('Original signal')
else
    minrz=min(min(tz));
    maxrz=max(max(tz));
    tzlab1=(minrz+min(T)/f-1)*f;
    tzlab2=(maxrz+min(T)/f-1)*f;
    y=c*s;
    yy=X0(minrz:maxrz,:);
    maxmax=1.1*max(max(y));
    
    subplot(1,2,1);
    plot(tzlab1:f:tzlab2,yy);
    axis([tzlab1 tzlab2 0 maxmax]);
    xlabel('elution time (s)')
    ylabel('intensity')
    title('Original signal')
    subplot(1,2,2);
    plot(tzlab1:f:tzlab2,y);
    axis([tzlab1 tzlab2 0 maxmax]);
    xlabel('elution time (s)')
    ylabel('intensity')
    title('Resolved signal')
    %plot(tzlab1:f:tzlab2,sum(y,2),'r','linewidth',2.0)
end
end
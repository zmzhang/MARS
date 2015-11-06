function mass_specstrum_plot(s,r,W)

%bar(W,s,'k')
bar(W,s,'r','barwidth',0.3);
[aa,index1]=sort(s,'descend');
xic1=index1(1:10)+min(W)-1;
text(xic1(1:10),aa(1:10),num2str([xic1(1:10)]'),'color','b')
text(800,0.95,['MSRT' num2str(r)]);
%set(gca,'tickdir','out')
xlabel('m/z')
ylabel('Intensity')
box on
% ax2 = axes('Position',get(gca,'Position'), 'XAxisLocation','top','YAxisLocation','right',...
%            'Color','none','XColor','k','YColor','k');
% set(ax2,'YTick',[]);
% set(ax2,'XTick',[]);
% box on
end


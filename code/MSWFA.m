function [region,s] = MSWFA(y,spec,max_pc,w,threshold)
%   SWFA:Subwindows Factor Analysis
%   Detailed explanation goes here

    L=size(y,1)-w+1;
    [~,~,F] = svd(spec);
    s=zeros(size(y,1),1);
    for j=1:L
       y1=y(j:j+w-1,:);
       [~,~,E] = svd(y1);
       M=E(:,1:max_pc)'*F(:,1)*F(:,1)'*E(:,1:max_pc);
       s(j+floor((w-1)/2))=max(svd(M));
    end
    vsel=find(s>=threshold);
    vse=(s>=threshold);
    if ~isempty(vsel);
        b=round(size(y,1)/2);
        [~,f]=min(abs(b-vsel));
        start=find(vse(1:vsel(f))==0, 1, 'last' );
        ends=vsel(f)+find(vse(vsel(f):size(y,1))==0, 1 )-1;
        region=start:ends;
    else
        region=[];
    end
end
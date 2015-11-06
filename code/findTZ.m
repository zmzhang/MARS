function [TZ,SEGMENTS] = findTZ(X0,S,R,w,threshold,pw)

[rows,~,nums]=size(X0);
[ns,~]=size(S);
max_pc=5;
SEGMENTS=zeros(nums,2,ns);
TZ=zeros(nums,2,ns);
half_pw=round(pw/2);

for i=1:ns
    tz=zeros(nums,2); % target zone
    segments=zeros(nums,2);
    for j=1:nums
        X=X0(:,:,j);
        s=S(i,:);
        rt=R(i);
        xblock=max([1 rt-half_pw-1]):min([rows rt+half_pw]);
        segments(j,:)=[xblock(1) xblock(end)];
        XX=X(xblock,:);
        %         total=sum(X(xblock,:));
        %         ve=find( total<0.001*max(total));
        %         XX(:,ve)=0;
        [region,eig] = MSWFA(XX,s,max_pc,w,threshold);
        if ~isempty(region)
            blocks=min(region):max(region);
            tzpoint=xblock(blocks);
            tz=[min(tzpoint) max(tzpoint)];
        end
    end
    
    SEGMENTS(:,:,i)=segments;
    TZ(:,:,i)=tz;
    disp([num2str(i),'/',num2str(ns),'th of findTZ ']);
end
end
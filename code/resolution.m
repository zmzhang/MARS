function  [C,AREAS,HEIGS,COEF] = resolution(X0,S,TZ,PCS,NEE)
[rows,cols,nums]=size(X0);
thre=0.95;
ns=size(TZ,3);
AREAS=zeros(nums,ns);
HEIGS=zeros(nums,ns);
C=cell(nums,ns);
COEF=cell(nums,ns);
for i=1:ns
    pw=TZ(:,:,i);
    pcs=PCS(:,i);
    nee=NEE(:,i);
    areas=zeros(nums,1);
    heigs=zeros(nums,1);
    s=S(i,:);
    coefs=[];
    c=[];
    for j=1:nums
        if  pw(j,1)==0
            COEF(j,i)={coefs};
            C(j,i)={c};
            areas(j)=0;
            heigs(j)=0;
        else
            needle=nee(j)-pw(j,1)+1;
            X=X0(:,:,j);
            y=X(pw(j,1):pw(j,2),:);
            [c,num] = ittfa(y,needle,pcs(j));
            %  polyfit
            for jr=1:size(y,2)
                R = corrcoef(c,y(:,jr)');
                r=R(1,2);
                coefs=[coefs abs(r)];
            end
%             [~,index]=sort(abs(coefs),'descend');
%             col=index(1:round(ratio*cols/100));
            col=find(coefs>=thre);
            if ~isempty(col)
            Z=c*s;
            k = polyfit(Z(:,col),y(:, col),1);
            yy = k(1)*Z;
            areas(j)=sum(sum(yy));
            heigs(j)=max(sum(yy,2));       
            cc=k(1)*c;
            else 
                areas(j)=0;
                heigs(j)=0;  
                cc=zeros(length(c),1);
            end
            %             figure;plot(sum(y,2)); hold on
            %             plot(sum(yy,2),'g');title(['conc profile',num2str(j)]);
            C(j,i)={cc};
            COEF(j,i)={coefs};
        end
        AREAS(:,i)=areas;
        HEIGS(:,i)=heigs;
        
        %close all
        disp([num2str(i),'/',num2str(ns),'th of resolution']);
    end
    
end
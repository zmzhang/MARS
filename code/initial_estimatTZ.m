function [PCS,NEE] = initial_estimatTZ(data,S,PW)
[rows,~,nums]=size(data);
ns=size(PW,3);
PCS=zeros(nums,ns);
NEE=zeros(nums,ns);
for i=1:ns
    pw=PW(:,:,i);
    needles=zeros(nums,1);
    pcs=zeros(nums,1);
    s=S(i,:);
    for j=1:nums
        
        if pw(j,1)~=0
            X=data(:,:,j);
            y=X(pw(j,1):pw(j,2),:);
            [A,D] = ferp_pcs(y,6); %   estimate the number of components of blocks
            %[A,D]=pc_number(y,10);
            %A = pc_number12(y,1.1);
            
            %  estimate the chromotograpgy profile (ittfa) and calculate the area
            %  confirm the initial needle position
            coefs=[];
            for it=1:size(y,1)
                coe=corrcoef(s,y(it,:));
                coefs=[coefs coe(1,2)];
            end
            [~,needle]=max(coefs);
            needles(j)=pw(j,1)+needle-1;
            pcs(j)=A;
        end
    end
    NEE(:,i)=needles;
    PCS(:,i)=pcs;
    disp([num2str(i),'/',num2str(ns),'th of initial_estimatTZ ']);
end
end
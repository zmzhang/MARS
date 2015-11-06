function RESULT = mars(X0,S,R,w,threshold,mpw)
%** MS-Asisted resolution of signal (MARS) - - - - *********************
%
%function manyareas = mars(X0,S,R,w,threshold,mpw);
%
%      where
%
% INPUT VALUES:
%
%       X0 : experimental data matrix (m x n x k)
%             m: number of chromatography; n: number of spectrum; k: number of samples
%         S : initial estimates of the spectrum profiles
%         R: initial position estimates of components
%         w: moving window (5 is the default)
%        pw: peak width (50 is default)
% threshold: eigenvalue (0.9 is the default)

% OUTPUT VALUES:
%
%  RESULT:  AREAS,HEIGS,TZ,PCS,NEE,segmengts,C;
%           AREAS: peak area
%           HEIGS: peak height
%           segmengts: small segments chromatographic profile located by Rt
%           TZ: target zone
%           PCS: components number of target zone
%           NEE: initial estimation of needle position
%           C: concentration profile
%********************.................................****.................*****

[TZ,SEGMENTS] = findTZ(X0,S,R,w,threshold,mpw);

[PCS,NEE] = initial_estimatTZ(X0,S,TZ);

[C,AREAS,HEIGS,COEF] = resolution(X0,S,TZ,PCS,NEE);

RESULT.AREAS=AREAS;RESULT.HEIGS=HEIGS;RESULT.PCS=PCS;
RESULT.NEE=NEE;RESULT.TZ=TZ;RESULT.SEGMENTS=SEGMENTS;RESULT.C=C;

end

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
            tz=[max([0 min(tzpoint)]) min([rows max(tzpoint)])];
        end
    end
    
    SEGMENTS(:,:,i)=segments;
    TZ(:,:,i)=tz;
    disp([num2str(i),'/',num2str(ns),'th of findTZ ']);
end
end

function [PCS,NEE] = initial_estimatTZ(data,S,TZ)
[rows,~,nums]=size(data);
ns=size(TZ,3);
PCS=zeros(nums,ns);
NEE=zeros(nums,ns);
for i=1:ns
    pw=TZ(:,:,i);
    needles=zeros(nums,1);
    pcs=zeros(nums,1);
    s=S(i,:);
    for j=1:nums
        
        if pw(j,1)~=0
            X=data(:,:,j);
            y=X(pw(j,1):pw(j,2),:);
            [A,D] = mars_pcs(y,6); %   estimate the number of components of blocks
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

function  [C,AREAS,HEIGS,COEF] = resolution(X0,S,TZ,PCS,NEE)
[rows,cols,nums]=size(X0);
thre=0.95;
ns=size(TZ,3);
AREAS=zeros(nums,ns);
HEIGS=zeros(nums,ns);
C=cell(nums,ns);
COEF=cell(nums,ns);
for i=1:ns
    tz=TZ(:,:,i);
    pcs=PCS(:,i);
    nee=NEE(:,i);
    areas=zeros(nums,1);
    heigs=zeros(nums,1);
    s=S(i,:);
    coefs=[];
    c=[];
    for j=1:nums
        if  tz(j,1)==0
            COEF(j,i)={coefs};
            C(j,i)={c};
            areas(j)=0;
            heigs(j)=0;
        else
            needle=nee(j)-tz(j,1)+1;
            X=X0(:,:,j);
            y=X(tz(j,1):tz(j,2),:);
            [c,num] = ittfa(y,needle,pcs(j));
            %  polyfit
            for jr=1:size(y,2)
                rr = corrcoef(c,y(:,jr)');
                r=rr(1,2);
                coefs=[coefs abs(r)];
            end
           %   [~,index]=sort(abs(coefs),'descend');
           %    col=index(1:round(ratio*cols/100));
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

end
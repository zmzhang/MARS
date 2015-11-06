function [X,T,W]=cdftomat(cdffilename,bin)
% matfilename=tempname;
% snc2mat(cdffilename,matfilename);
% eval(['load ' matfilename]);

filename=cdffilename;
T=ncread(filename,'scan_acquisition_time');
TIC=ncread(filename,'total_intensity');
scan_index=ncread(filename,'scan_index');
mass_range_min=ncread(filename,'mass_range_min');
mass_range_max=ncread(filename,'mass_range_max');
mass_values=ncread(filename,'mass_values');
intensity_values=ncread(filename,'intensity_values');
resolution=ncread(filename,'resolution');

% minibin=2*min(resolution(1));%it sees bin size must than 3*
% disp(['The m/z resolution is ' num2str(resolution.data(1))]);
% if bin<minibin
%     warning(['The minimum bin size is ' num2str(minibin)]);
%     return;
% end


N=length(scan_index)-1;
T=T(1:N);
masslow=min(mass_range_min);
masshigh=max(mass_range_max);
W=[masslow:bin:masshigh]';
M=length(W);
X=zeros(N,M);
%figure(1);
for i=1:N
    msnext=scan_index(i+1);
    if msnext>length(mass_values)
        T(i:N)=[];
        X(i:N,:)=[];
        break;
    else
        ptind=scan_index(i)+1:scan_index(i+1);
        ms=mass_values(ptind);
        sp=intensity_values(ptind);
        %         plot(ms,sp);
        %         title(T(i));
        %         pause
        %ind=floor((ms-masslow)/bin)+1;%the column index %one can not use round here
        ind=round((ms-masslow)/bin)+1;
        position= ind>M;
        ind(position)=[];
        sp(position)=[];
        [ind2]=unique(ind);%unique index
        if length(ind2)~=length(ind)
            sp2=zeros(1,length(ind2));%new variable with the same size of new index
            for j=1:length(ind2)
                tempind=find(ind==ind2(j));
                sp2(j)=sum(sp(tempind));
            end
        else
            sp2=sp;
            ind2=ind;
        end
        X(i,ind2)=sp2;
    end
end



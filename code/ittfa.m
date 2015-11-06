function [cout,num] = ittfa(y,needle,A,kesse)
   if nargin<4;kesse=10^-6;end
   [u,s,~]=svd(y);
   T=u*s;
   T=T(:,1:A);
   row=size(y,1);
   cin=zeros(row,1);
   cin(needle)=1;
   cout=cin;

   for iter=1:1000
      vect=cout;
      cout=T*pinv(T'*T)*T'*cout;
      cout(cout<0)=0;
      [cout]=unimod(cout,1.1,2); 
      cout=cout/norm(cout);
      %cout(1)=0;
      
      kes=norm(cout-vect);  
      if kes<kesse||iter==1000
         num=iter;
         break;
      end
   end  
   
end
  
  
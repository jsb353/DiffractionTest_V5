function TopDownMerge(A,iBegin, iMiddle,iEnd, B,column)

i=iBegin; j=iMiddle;

for k=iBegin:iEnd
    
   if (i<iMiddle&&(j>iEnd || A(i,column)>A(j,column)))
      B(k,:)=A(i,:);
      i=i+1;       
   else
       B(k,:)=A(j,:);
       j=j+1;       
   end    
end

end
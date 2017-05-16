function TopDownSplitMerge(B,iBegin,iEnd, A,column)

if (iEnd - iBegin<2)
    
else
    iMiddle=floor((iEnd+iBegin)/2);
    TopDownSplitMerge(A,iBegin, iMiddle,B,column);
    TopDownSplitMerge(A, iMiddle,iEnd, B,column);
    TopDownMerge(B,iBegin,iMiddle,iEnd,A,column);
    
end

end
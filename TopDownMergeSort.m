function B=TopDownMergeSort(A,B,n,column)

B=A;
TopDownSplitMerge(B,1,n,A,column);

end
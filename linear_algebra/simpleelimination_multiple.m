function solution=simpleelimination_multiple(A,k)
%EP501 Homework 1
%Ex. 1
%Part a

%This function provides a simple forward elimination method as already
%implemented in class examples that can be used with any matrix A and
%multiple RHS. k is a right hand side values vector(columns)
[o,j]=size(k);
[r,p]=size(A);
for i=1:j
    b=k(:,i);
    nref=length(b);    %system size for reference problem
    f=nref*(i-1);
    y=p*(i-1);
    %note that the elimination procedure coded below modifies the matrix B
    Awork=cat(2,A,b);   %This is our working version of the matrix used to perform elimination (i.e. it will be modified)
    solution(1+f,:)=Awork(1,:);
    for ir1=2:nref       %loop over rows from 2 to n performing elimination, this index marks what row we are starting the elimination from (i.e. using) for this particular column
        for ir2=ir1:nref   %this index marks the present position where elimination is being performed - i.e. where we are applying the elementary row operations
            fact=Awork(ir2,ir1-1);     %multiplier of the variable we are attempting to eliminate, its ir-1 column of this row
            Awork(ir2,:)=Awork(ir2,:)-fact/Awork(ir1-1,ir1-1).*Awork(ir1-1,:);    %subtract off previous row modified by a factor that eliminates the ir-1 column term in this row (so it has only super-diagonal elements), this is a little bit wasteful as it uses entire row...
            solution(ir2+f,:)=Awork(ir2,:);
        end %for
    end %for
end


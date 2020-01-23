function D=mahal_sam(X,M)

n=size(X,1);
D=zeros(n);
for i=1:n
    for j=i+1:n
        D(i,j)=X(i,:)*M*X(j,:)';
    end
end
D=D+D';

end
function P = paralInv(Sigma)

clc
Sigma=rand(8,8);
p=size(Sigma,1);
nc=3;
nr=3;
pc=floor(p/nc)+1;
pr=floor(p/nr)+1;

row1=1;
for i=1:nr
    if row1+pr-1<p
        rows=[row1:row1+pr-1];
    else
        rows=[row1:p];
    end
    display(rows);
    C{i}=Sigma(rows,:);
    row1=row1+pr;
end

col1=1;
for j=1:nc
    if col1+pc-1<p
        cols=[col1:col1+pc-1];
    else
        cols=[col1:p];
    end
    display(cols);
    Et=eye(p);
    Et=Et(:,cols);
    row1=1;
    for i=1:nr
        if row1+pr-1<p
            rows=[row1:row1+pr-1];
        else
            rows=[row1:p];
        end
        display(rows);
        E{i}=Et(rows,:);
        row1=row1+pr;
    end
    
    col1=col1+pc;
end



end
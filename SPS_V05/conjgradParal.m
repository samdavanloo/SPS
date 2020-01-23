function [x] = conjgradParal(A,b,x,nr)
    r=b-A*x;
    p=r;
    rsold=r'*r;
    
    d=size(A,1);
    if floor(d/nr)<(d/nr)
        dr=floor(d/nr)+1;
    else
        dr=d/nr;
    end

    row1=1;
    for i=1:nr
        if row1+dr-1<d
            rows=[row1:row1+dr-1];
        else
            rows=[row1:d];
        end
        Ar{i}=A(rows,:);
        row1=row1+dr;
    end
 
    for i=1:1e6
        parfor j=1:nr %parfor
            Apr{j}=Ar{j}*p;
        end
        Ap=[];
        for j=1:nr
            Ap=[Ap;Apr{j}];
        end
        alpha=rsold/(p'*Ap);
        x=x+alpha*p;
        r=r-alpha*Ap;
        rsnew=r'*r;
        if sqrt(rsnew)<1e-10
              break;
        end
        p=r+rsnew/rsold*p;
        rsold=rsnew;
    end
end
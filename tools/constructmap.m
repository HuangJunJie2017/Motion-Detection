function map = constructmap(x,y)
    map =zeros(x*y,2);
    for i =1:x
        for j =1:y
            idx = (i-1)*y+j;
            map(idx,1)=i;
            map(idx,2)=j;
        end
    end
end
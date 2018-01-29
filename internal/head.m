function res=head(list)

% shows the first 6 lines

if size((list),1)>1 & size((list),1)>6 & size((list),2)>1
    for i=1:6
        res(i,:)=list(i,:);
    end
end

if size((list),1)>1 & size((list),1)<6 & size((list),2)>1
   
    for i=1:size((list),1)
        res(i,:)=list(i,:);
    end
end

if size((list),1)>1 & size((list),1)<6 & size((list),2)==1
   
    for i=1:size((list),1)
        res(i,:)=list(i,:);
    end
end


if size((list),1)>1 & size((list),1)>=6 & size((list),2)==1
   
    for i=1:6
        res(i,:)=list(i,:);
    end
end


if size((list),1)==1 & size((list),2)>1 & size((list),2)<6
   
    for i=1:size((list),2)
        res(i,:)=list(i,:);
    end
end

if size((list),1)==1 & size((list),2)>1 & size((list),2)>=6
   
    for i=1:6
        res(i,:)=list(:,i);
    end
end

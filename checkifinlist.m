function checkedlist = checkifinlist(sublist,superlist)

checkedlist=[]; k = size(sublist,2);
for i = 1:size(sublist,1)
    if prod(ismember(nchoosek(sublist(i,:),k-1),superlist,'rows'))==1
        checkedlist = [checkedlist; sublist(i,:)];
    end
end

function [treeIndex,Dist] = tree_delete_repeats(treeValues,searchValues)
%% [treeIndex,Dist] = tree_delete_repeats(treeValues,searchValues)
%treeIndex(s) is the index of the closest value in treeValues to
%searchValues(s) such that each value in treeValues can be a neighbour to
%atmost 1 value of searchValues.
%Dist(s) is the corresponding distance.

tree = KDTreeSearcher(treeValues);
[treeIndex,Dist] = tree.knnsearch(searchValues);
RepeatNearest_treeIndex = 1;
while ~isempty(RepeatNearest_treeIndex)
    treeIndexCount = count(treeIndex);
    RepeatNearest_treeIndex = unique(treeIndexCount(treeIndexCount(:,2)>=2,1))';
    treeIndicesNotUsed = setdiff(1:size(treeValues,1),treeIndex);
    treeNotUsed = KDTreeSearcher(treeValues(treeIndicesNotUsed,:));
    for s=RepeatNearest_treeIndex
        Set = find(treeIndex==s);
        [MinDist,MinIndex] = min(Dist(Set));
        %Change neighbour for the non closest neighbours
        if sum(Dist(Set)==MinDist)>1
            %If two have same distance, change first one
            Dist(Set(MinIndex))=inf;
        end
        searchValues_change_index = Set(Dist(Set)>MinDist);
        [IndexNewTree,Dist(searchValues_change_index)] = ...
            treeNotUsed.knnsearch(searchValues(searchValues_change_index,:));
        treeIndex(searchValues_change_index) = treeIndicesNotUsed(IndexNewTree);
    end
end
end


function C=count(A)
%C(:,1)=A
%C(s,2) is equal to 1 if it is the first time that element has appeared, 2
%if second, 3 if third etc. 
A=A(:);
[C,i]=sort(A);
j=find([true;diff(C)]);
lgt=diff(j);
c=ones(size(A));
c(j(2:end))=1-lgt;
c=cumsum(c);
c(i)=c;
C=[A,c];
end


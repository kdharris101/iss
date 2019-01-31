function iss_show_genes(GenesToShow, NotThese, Text)
% ShowGenes(GenesToShow, NotThese, Text)
% 
% with in situ sequencing plots, show genes with names in GenesToShow
% but not also in NotThese. 
% 
% If GenesToShow is empty, defaults to all
% If NotThese is empty, or not provided, defaults to none
% If Text is 'off', don't show celltype labels (default is 'on')

if nargin<3
    Text='on';
end

chdn0 = get(gca, 'children');
MyChildren = false(size(chdn0));
MyText = MyChildren;
for i=1:length(chdn0)
    if strcmp(chdn0(i).Type, 'line') | strcmp(chdn0(i).Type, 'scatter')
        %MyChildren = [MyChildren; i];
        MyChildren(i)=1;
    elseif strcmp(chdn0(i).Type, 'text')
        MyText(i)=1;
    end
end


set(chdn0(MyText), 'Visible', Text);


chdn = chdn0(MyChildren); % the last one is the image
GeneNames = get(chdn, 'DisplayName');

if nargin<1 | isempty(GenesToShow) 
    GenesToShow = GeneNames;
end

if nargin>=2
    GenesToShow = setdiff(GenesToShow, NotThese);
end


vIDs = find(ismember(GeneNames, GenesToShow));
nvIDs = setdiff(1:length(chdn),vIDs);

set(chdn(vIDs), 'Visible', 'on');
set(chdn(nvIDs), 'Visible', 'off');


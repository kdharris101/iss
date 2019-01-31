function c = ca1_pie_colors(o)
% creates a cell array to store in o.ClassCollapse to nicely display the
% colors of CA1 cells

% how much do individual classes vary compared to the group average
NoiseSize = .2;

randn('state', 1);

% copied from change_gene_symbols
non_neuron = hsv2rgb([0 0 1]);
pc_or_in =   hsv2rgb([.4 .5 .5]);
less_active =   hsv2rgb([.3 .2 .7]);
pc =        hsv2rgb([1/3 1 1]);
pc2 =       hsv2rgb([.27 1 .7]);
in_general = hsv2rgb([2/3 1 1]);

sst =   hsv2rgb([.55 1 1]);
pvalb = hsv2rgb([.7 .8 1]);
ngf =   hsv2rgb([.85 1 1]);
cnr1 =  hsv2rgb([ 1 1 1]);
vip =   hsv2rgb([ .13 1 1]);
cxcl14= hsv2rgb([.1 1 .6]);

ivy  =   hsv2rgb([.85 .5 .6]);

% dictionary of colors
d = {{'Sst'}, sst ; {'Pvalb'}, pvalb ; {'Cacna2d1.Lhx6'}, ivy ; ...
    {'Cacna2d1.Ndnf'}, ngf ; {'Ntng1'}, pc_or_in; ...
    {'Cck.Cxcl14'}, cxcl14 ; {'Cck.Lmo1', 'Cck.Calca', 'Cck.Sema5a'}, cnr1; ...
    {'Calb2', 'Vip'}, vip};

% output
c = cell(0,3);

for i=1:size(d,1)
    MyClasses = [];
    for j=1:length(d{i,1})
        MyClasses = union(MyClasses, strmatch(d{i,1}{j}, o.ClassNames));
    end
    for j=1:length(MyClasses)
        cn = o.ClassNames(MyClasses(j));
        color = d{i,2}*(1-NoiseSize) + NoiseSize.*rand(1,3);
        color(color<0) = 0; color(color>1) = 1;
        c = vertcat(c, {cn, cn{1}, color});
    end
end
    
% manually-defined classes 
c = vertcat(c, {{'PC.CA1'}, 'PC CA1', pc});
c = vertcat(c, {{'PC.CA2', 'PC.CA3'}, 'PC other', less_active});
c = vertcat(c, {{'Astro', 'Endo', 'Oligo', 'Eryth', 'Vsmc', 'Microglia', 'Choroid'}...
    , 'Non neuron', [.5 .5 .5]});
c = vertcat(c, {{'Zero'}, 'Uncalled', [.1 .1 .1]});

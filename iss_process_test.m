o = iss;

o.FileBase = cell(6,1);
o.FileBase{1} = '170315_161220KI_4-1_b1';
o.FileBase{2} = '170315_161220KI_4-1_b2';
o.FileBase{3} = '170315_161220KI_4-1_b3';
o.FileBase{4} = '170315_161220KI_4-1_b4';
o.FileBase{5} = '170315_161220KI_4-1_b5';
o.FileBase{6} = '170315_161220KI_4-1_SW';

o = o.extract_and_filter;

if 0
    o.EmptyTiles = ones(size(o.EmptyTiles));
    o.EmptyTiles(8:9,5:6)=0;
end

save oExtract o -v7.3;

o = o.register;
save oRegister o -v7.3;

o = o.find_spots;
save oFind_spots o -v7.3;

o = o.call_spots;
save oCall_spots o -v7.3;

o = o.add_extra_genes;
save oExtra_genes o -v7.3;

figure(100);
o.plot();

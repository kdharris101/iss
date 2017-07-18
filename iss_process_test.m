o = iss;

o.FileBase = cell(6,1);
o.FileBase{1} = '170315_161220KI_4-1_b1';
o.FileBase{2} = '170315_161220KI_4-1_b2';
o.FileBase{3} = '170315_161220KI_4-1_b3';
o.FileBase{4} = '170315_161220KI_4-1_b4';
o.FileBase{5} = '170315_161220KI_4-1_b5';
o.FileBase{6} = '170315_161220KI_4-1_SW';

o = o.extract_and_filter;
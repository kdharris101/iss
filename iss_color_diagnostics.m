function iss_color_diagnostics(o)
% iss_color_diagnostics(o)
%
% does some diagnostics of signal and noise in each color channel and
% round. Makes a matrix of histograms. In each one, the blue bars are the
% histogram of intensities on spots which should not have fluorescence on
% that color and round. The red bars are the ones that should have
% fluorescence. So the red bars should be to the right of the blue ones.
% Note the y-axis is on a log scale .

nc = size(o.cSpotColors,2);
nr = size(o.cSpotColors,3);

figure(349075); clf; 

for c=1:nc
    for r=1:nr

        rc = (r-1)*o.nBP + c;
        subplot(nr,nc,rc);


        ShouldBe1 = (o.UnbledCodes(o.SpotCodeNo,rc)>0);

        cla; hold on
        histogram(o.cSpotColors(~ShouldBe1,c,r), 'FaceColor', 'b', 'EdgeColor', 'b');
        histogram(o.cSpotColors(ShouldBe1,c,r), 'FaceColor', 'r', 'EdgeColor', 'r');
        
        title(sprintf('Color %d Round %d', c-1, r))
        set(gca, 'yscale', 'log')
    end
end
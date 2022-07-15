function I = colormapTestImage(map)
% colormapTestImage Create or display colormap test image.
%   I = colormapTestImage creates a grayscale image matrix that is useful
%   for evaluating the effectiveness of colormaps for visualizing
%   sequential data. In particular, the small-amplitude sinusoid pattern at
%   the top of the image is useful for evaluating the perceptual uniformity
%   of a colormap.
%
%   colormapTestImage(map) displays the test image using the specified
%   colormap. The colormap can be specified as the name of a colormap
%   function (such as 'parula' or 'jet'), a function handle to a colormap
%   function (such as @parula or @jet), or a P-by-3 colormap matrix.
%
%   EXAMPLES
%
%     Compute the colormap test image and save it to a file.
%
%       mk = colormapTestImage;
%       imwrite(mk,'test-image.png');
%
%     Compare the perceptual characteristics of the parula and jet
%     colormaps.
%
%       colormapTestImage('parula')
%       colormapTestImage('jet')
%
%   NOTES
%
%   The image is inspired by and adapted from the test image proposed in
%   Peter Kovesi, "Good Colour Maps: How to Design Them," CoRR, 2015,
%   https://arxiv.org/abs/1509.03700
%
%   The upper portion of the image is a linear ramp (from 0.05 to 0.95)
%   with a superimposed sinusoid. The amplitude of the sinusoid ranges from
%   0.05 at the top of the image to 0 at the bottom of the upper portion.
%
%   The lower portion of the image is a pure linear ramp from 0.0 to 1.0.
%
%   This test image differs from Kovesi's in three ways:
%
%     (a) The Kovesi test image superimposes a sinusoid on top of a
%     full-range linear ramp (0 to 1). It then rescales each row
%     independently to have full range, resulting in a linear trend slope
%     that slowly varies from row to row. The modified test image uses the
%     same linear ramp (0.05 to 0.95) on each row, with no need for
%     rescaling.
%
%     (b) The Kovesi test image has exactly 64 sinusoidal cycles
%     horizontally. This test image has 64.5 cycles plus one sample. With
%     this modification, the sinusoid is at the cycle minimum at the left
%     of the image, and it is at the cycle maximum at the right of the
%     image. With this modification, the top row of the modified test image
%     varies from exactly 0.0 on the left to exactly 1.0 on the right,
%     without rescaling.
%
%     (c) The modified test image adds to the bottom of the image a set of
%     rows containing a full-range (0.0 to 1.0) linear ramp with no
%     sinusoidal variation. That makes it easy to view how the colormap
%     appears with a full-range linear ramp.
%
%   Reference: Peter Kovesi, "Good Colour Maps: How to Design Them,"
%   CoRR, 2015, https://arxiv.org/abs/1509.03700
%   Steve Eddins
%   Copyright 2017 The MathWorks, Inc.
% Compare with 64 in Kovesi 2015. Adding a half-cycle here so that the ramp
% + sinusoid will be at the lowest part of the cycle on the left side of
% the image and at the highest part of the cycle on the right side of the
% image.
num_cycles = 64.5;
if nargin < 1
    I = testImage(num_cycles);
else
    displayTestImage(map,num_cycles)
end
function I = testImage(num_cycles)
pixels_per_cycle = 8;
A = 0.05;
% Compare with width = pixels_per_cycle * num_cycles in Kovesi 2015. Here,
% the extra sample is added to fully reach the peak of the sinusoid on the
% right side of the image.
width = pixels_per_cycle * num_cycles + 1;
% Determined by inspection of 
% http://peterkovesi.com/projects/colourmaps/colourmaptest.tif
height = round((width - 1) / 4);
% The strategy for superimposing a varying-amplitude sinusoid on top of a
% ramp is somewhat different from Kovesi 2015. For each row of the test
% image, Kovesi adds the sinusoid to a full-range ramp and then rescales
% the row so that ramp+sinusoid is full range. A benefit of this approach
% is that each row is full range. A drawback is that the linear trend of
% each row varies as the amplitude of the superimposed sinusoid changes.
%
% Our approach here is a modification. The same linear ramp is used for
% every row of the test image, and it goes from A to 1-A, where A is the
% amplitude of the sinusoid. That way, the linear trend is identical on
% each row. The drawback is that the bottom of the test image goes from
% 0.05 to 0.95 (assuming A = 0.05) instead of from 0.00 to 1.00.
ramp = linspace(A, 1-A, width);
k = 0:(width-1);
x = -A*cos((2*pi/pixels_per_cycle) * k);
% Amplitude of the superimposed sinusoid varies with the square of the
% distance from the bottom of the image.
q = 0:(height-1);
y = ((height - q) / (height - 1)).^2;
I1 = (y') .* x;
% Add the sinusoid to the ramp.
I = I1 + ramp;
% Add region to the bottom of the image that is a full-range linear ramp.
I = [I ; repmat(linspace(0,1,width), round(height/4), 1)];
function displayTestImage(map,num_cycles)
name = '';
if isstring(map)
    map = char(map);
    name = map;
    f = str2func(map);
    map = f(256);
elseif ischar(map)
    name = map;
    f = str2func(map);
    map = f(256);
elseif isa(map,'function_handle')
    name = func2str(map);
    map = map(256);
end
I = testImage(num_cycles);
[M,N] = size(I);
% Display the image with a width of 2mm per cycle.
display_width_cm = num_cycles * 2 / 10;
display_height_cm = display_width_cm * M / N;
fig = figure('Visible','off',...
    'Color','k');
fig.Units = 'centimeters';
% Figure width and height will be image width and height plus 2 cm all the
% way around.
margin = 2;
fig_width = display_width_cm + 2*margin;
fig_height = display_height_cm + 2*margin;
fig.Position(3:4) = [fig_width fig_height];
ax = axes('Parent',fig,...
    'DataAspectRatio',[1 1 1],...
    'YDir','reverse',...
    'CLim',[0 1],...
    'XLim',[0.5 N+0.5],...
    'YLim',[0.5 M+0.5]);
ax.Units = 'centimeters';
ax.Position = [margin margin display_width_cm display_height_cm];
ax.Units = 'normalized';
box(ax,'off')
im = image('Parent',ax,...
    'CData',I,...
    'XData',[1 N],...
    'YData',[1 M],...
    'CDataMapping','scaled');
if ~isempty(name)
    title(ax,name,'Color',[0.8 0.8 0.8],'Interpreter','none')
end
% Draw scale line.
pixels_per_centimeter = N / (display_width_cm);
x = [0.5 5*pixels_per_centimeter];
y = (M + 30) * [1 1];
line('Parent',ax,...
    'XData',x,...
    'YData',y,...
    'Color',[0.8 0.8 0.8],...
    'Clipping','off');
text(ax,mean(x),y(1),'5cm',...
    'VerticalAlignment','top',...
    'HorizontalAlignment','center',...
    'Color',[0.8 0.8 0.8]);
colormap(fig,map)
fig.Visible = 'on';

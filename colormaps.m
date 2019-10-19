function  x=colormaps(u)
% x = colormaps
% return a list of builting colour maps, as a cell string 
% ( 18 items )
%
% x = colormaps( i )
% returns the actual colourmap with index i
% as a matrix of n x 3 colour values.

 x={'hsv' %        - Hue-saturation-value color map.
    'hot' %        - Black-red-yellow-white color map.
    'gray' %       - Linear gray-scale color map.
    'bone' %       - Gray-scale with tinge of blue color map.
    'copper' %     - Linear copper-tone color map.
    'pink' %       - Pastel shades of pink color map.
    'white' %      - All white color map.
    'flag' %       - Alternating red, white, blue, and black color map.
    'lines' %      - Color map with the line colors.
    'colorcube' %  - Enhanced color-cube color map.
    'vga' %        - Windows colormap for 16 colors.
    'jet' %        - Variant of HSV.
    'prism' %      - Prism color map.
    'cool' %       - Shades of cyan and magenta color map.
    'autumn' %     - Shades of red and yellow color map.
    'spring' %     - Shades of magenta and yellow color map.
    'winter' %     - Shades of blue and green color map.
    'summer' %     - Shades of green and yellow color map.
}; 
if exist('u','var')
  x=x{u};
  x=colormap(x);
end
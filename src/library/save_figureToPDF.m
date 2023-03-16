
function file_baseName = save_figureToPDF( handle, file_baseName, fig_area, paper_DPI )

% Usage: save_figureToPDF( handle, fileBaseName, fig_area, paper_DPI)
%
% Saves a given figure specified by 'handle' into PDF format file with the name 'fileBaseName'.
% You can also specify the area of figure which is converted to a PDF file.
% If you want to adjust the resolution, it is possible to set the pixes per inch in PDF file.
%
% handle:        Handle to a Matlab figure
%
% fileBaseName:  The name of PDF file, not including '.pdf' extension.
%
% 'fig_area':    (optional) The area to be saved into PDF file.
%                Its format is [X(left_bottom) Y(left_bottom) X(right_top) Y(right_top)]
%                The unit is a 'normalized' unit which ranges from 0 to 1.
%                If this parameter isn't specified,
%                the whole area of figure will be printed to PDF file.
%
% 'paper_DPI':   (optional) Pixels per inch in PDF file.
%                If this value isn't specified, DPI of screen will be used instead.
%
% EXAMPLE:
%    save_figureToPDF( handle, 'figures/plot1' );
%    save_figureToPDF( handle, 'fig_1', [.2 .1 .8 .9] );
%    save_figureToPDF( handle, 'fig_2', [.2 .1 .8 .9], 200 );
%
% Copyright (C) Daeseob Lim (daeseob at gmail dot com), 2013.

	if ( ~exist( 'fig_area', 'var' ) || isempty(fig_area) )
		fig_area = [0 0 1 1];
	end
	assert( fig_area(1) <= fig_area(3) && fig_area(2) <= fig_area(4) );

	pixel_pos = getpixelposition( handle );
	% NOTE) If we cut the figure according to 'fig_area',
	%	the right and the bottom portions of the figure are fully shown.
	%	So, we give a little margin to right and bottom areas.
	deltaX_perPixel	= 1 / pixel_pos(3);
	deltaY_perPixel	= 1 / pixel_pos(4);
	hor_margin	= 2;	% in pixels
	vert_margin	= 2;	% in pixels

	h_pos	= fig_area(1) - hor_margin * deltaX_perPixel / 2;
	v_pos	= fig_area(2) - vert_margin * deltaY_perPixel / 2;
	width	= fig_area(3) - fig_area(1) + hor_margin * deltaX_perPixel;
	height	= fig_area(4) - fig_area(2) + vert_margin * deltaY_perPixel;

	% Determine the size of output image in terms of paper size (in cm)
	screen_DPI = get( 0, 'ScreenPixelsPerInch' );
	wholeFigureSize	= pixel_pos(3:4) / screen_DPI * 2.54;					% in cm
	validFigureSize	= wholeFigureSize .* [width height];	% in cm

	paperSize		= validFigureSize;
	leftMargin		= wholeFigureSize(1) * h_pos;
	bottomMargin	= wholeFigureSize(2) * v_pos;
	paperPosition	= [-leftMargin, -bottomMargin, wholeFigureSize];

	% Set the proper properties
	set( handle, 'PaperUnits', 'centimeters' );
	set( handle, 'PaperSize', paperSize );
	set( handle, 'PaperPosition', paperPosition );
	set( handle, 'InvertHardcopy', 'off' );
	
	% Print out to a PDF file (NOTE: Without '-r' option, default DPI is set to 600.)
	%full_name_EPS = [file_baseName '.eps'];
	full_name_PDF = [file_baseName '.pdf'];

	% NOTE) Edited by Danny on Nov 2, 2015.
	% There are two options for exporting a figure to PDF file.
	%	1) Calling print() function embedded in Matlab.
	%		. This is a default option. But, since Matlab 2014b, the new graphic system in Matlab
	%			has failed to export a high-quality PDF file just as seen on the screen.
	%			One critical issue was that it didn't print thin lines (<= 0.5 line width).
	%		. But, this problem seems to be solved in Matlab 2015b.
	%			So, I decided to come back to the use of print() function.
	%	2) Using print2eps() function provided by third-party library
	%		. This 'export_fig' library is so popular in Matlab community.
	%			It is known for high quality PDF export functionality.
	%			But, it seems it still suffers from failure in printing thin lines.
	%			Even worse, it fails to print 'imagesc()' map when the thickness of row vectors is not even.
	if ( exist('paper_DPI','var') )
		print( handle, full_name_PDF, '-dpdf', ['-r' num2str(paper_DPI)] );
	else
		default_DPI = 600;
		print( handle, full_name_PDF, '-dpdf', ['-r' num2str(default_DPI)] );
	end
	%{
	print2eps( full_name_EPS, handle );
	eps2pdf( full_name_EPS, full_name_PDF );
	delete( full_name_EPS );
	%}

end

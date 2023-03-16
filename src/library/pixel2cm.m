
function values_in_cm = pixel2cm( values_in_pixel )
	screen_DPI		= get( 0, 'ScreenPixelsPerInch' );
	cm_per_inch		= 2.54;
	values_in_cm	= values_in_pixel / screen_DPI * cm_per_inch;
end

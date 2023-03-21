
function values_in_pixel = cm2pixel( values_in_cm )
	screen_DPI		= get( 0, 'ScreenPixelsPerInch' );
	cm_per_inch		= 2.54;
	values_in_pixel	= values_in_cm / cm_per_inch * screen_DPI;
end

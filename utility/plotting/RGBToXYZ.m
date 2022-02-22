
function[xyz] = RGBToXYZ(rgb)

  r = rgb(1); g = rgb(2); b = rgb(3);

  % The following performs a "gamma correction" specified by the sRGB color
  % space.  sRGB is defined by a canonical definition of a display monitor and
  % has been standardized by the International Electrotechnical Commission (IEC
  % 61966-2-1).  The nonlinearity of the correction is designed to make the
  % colors more perceptually uniform.  This color space has been adopted by
  % several applications including Adobe Photoshop and Microsoft Windows color
  % management.  OpenGL is agnostic on its RGB color space, but it is reasonable
  % to assume it is close to this one.
  if ( r > 0.04045 ), r = (( r + 0.055 ) / 1.055)^2.4;
  else                r = r / 12.92; end
  if ( g > 0.04045 ), g = (( g + 0.055 ) / 1.055)^2.4;
  else                g = g / 12.92; end
  if ( b > 0.04045 ), b = (( b + 0.055 ) / 1.055)^2.4;
  else                b = b / 12.92; end

  %Observer. = 2 deg, Illuminant = D65
  x = r * 0.4124 + g * 0.3576 + b * 0.1805;
  y = r * 0.2126 + g * 0.7152 + b * 0.0722;
  z = r * 0.0193 + g * 0.1192 + b * 0.9505;
  
  xyz = [x y z];
end

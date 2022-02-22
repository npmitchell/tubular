function[rgb] = XYZToRGB(xyz)
  
  %ref_X = 0.9505;        %Observer = 2 deg Illuminant = D65
  %ref_Y = 1.000;
  %ref_Z = 1.089;
  
  x = xyz(1); y = xyz(2); z = xyz(3);
  r = x *  3.2406 + y * -1.5372 + z * -0.4986;
  g = x * -0.9689 + y *  1.8758 + z *  0.0415;
  b = x *  0.0557 + y * -0.2040 + z *  1.0570;

  % The following performs a "gamma correction" specified by the sRGB color
  % space.  sRGB is defined by a canonical definition of a display monitor and
  % has been standardized by the International Electrotechnical Commission (IEC
  % 61966-2-1).  The nonlinearity of the correction is designed to make the
  % colors more perceptually uniform.  This color space has been adopted by
  % several applications including Adobe Photoshop and Microsoft Windows color
  % management.  OpenGL is agnostic on its RGB color space, but it is reasonable
  % to assume it is close to this one.
  if (r > 0.0031308), r = 1.055 * r^( 1 / 2.4 ) - 0.055;
  else r = 12.92 * (r); end
  if (g > 0.0031308), g = 1.055 * g^( 1 / 2.4 ) - 0.055;
  else  g = 12.92 * (g); end
  if (b > 0.0031308), b = 1.055 * b^( 1 / 2.4 ) - 0.055;
  else b = 12.92 * (b); end

  % Clip colors. ideally we would do something that is perceptually closest
  % (since we can see colors outside of the display gamut), but this seems to
  % work well enough.
  maxVal = r;
  if (maxVal < g), maxVal = g; end
  if (maxVal < b), maxVal = b; end
  if (maxVal > 1.0)    
    r = r/maxVal;
    g = g/maxVal;
    b = b/maxVal;
  end
  if (r<0), r=0; end
  if (g<0), g=0; end
  if (b<0), b=0; end
  
  rgb = [r g b];
end

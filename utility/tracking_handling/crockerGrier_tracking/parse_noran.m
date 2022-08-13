function [nx,ny,nz,xdimension,ydimension,zdimension] = parse_noran(filename)
  
  
  fid = fopen(filename);
  head = char(fread(fid,3200))';
  
%  protocol     = get_item(head,'Protocol',/string)
  
  objlensmag   = get_item(head,'ObjLensMag',0);
  objlensna    = get_item(head,'ObjLensNA',0);

  xlenscalibrt = get_item(head,'XLensCalibrt',0);
  ylenscalibrt = get_item(head,'YLensCalibrt',0);
  
  xdimension   = get_item(head,'Xdimension',0);
  ydimension   = get_item(head,'Ydimension',0);
  zdimension   = get_item(head,'Zdimension',0);

  irisvalue    = get_item(head,'IrisValue',0);
  slitwidth    = get_item(head,'SlitWidth',0);
  
  nx           = get_item(head,'WIDTH',0);
  ny           = get_item(head,'HEIGHT',0);
  nz           = get_item(head,'DIR_COUNT',0);
  
  zoomvalue    = get_item(head,'ZoomValue',0);
  adcdwell     = get_item(head,'ADCDwell',0);
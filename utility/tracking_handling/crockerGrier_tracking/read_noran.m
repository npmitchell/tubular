function a = read_noran(filename,numz)
  
  [nx,ny,nz,xdim,ydim,zdim] = parse_noran(filename);
  if nz == 0
      nz = 1;
  end
  
  size = nx * ny * nz;
  fid = fopen(filename);
  [filename,permission,machineformat] = fopen(fid);
  
  fseek(fid,0,'eof');
  f_size = ftell(fid);
  sz = f_size-size;
  frewind(fid);
  
  fseek(fid,sz,0);
  ftell(fid);
  if nargin == 2
    nz = numz;
    mat = repmat(uint8(0),[nx ny nz]);
  else    
    %nz = 0;
    mat = zeros(nx,ny);
  end
  
%frewind(fid);
%mat = multibandread(filename,[nx,ny,nz],'double',sz,'bip',machineformat);
if nz > 0
    for k=1:nz
%        mat(:,:,:) = fread(fid,[nx,ny,nz],'uchar');
        mat(:,:,k) = fread(fid,[nx,ny],'uchar');
        end
    else
    mat = fread(fid,[nx,ny],'uchar');
end

  fclose(fid);
  strcat(num2str(xdim),'  x  ',num2str(ydim),'  x  ',num2str(zdim),' um')
  a = mat;
    
function a = read_gdf(filename)
  
  % guess first that it is not natively written on a mac
  fid = fopen(filename,'r','l');
  if fid == -1 
    'not a real file'
    return
  end
  % first check the endian (I think this is correct)
  a = fread(fid,1,'ushort');
  stat = fseek(fid,0,-1);
  if a(1) == 256  % written on the mac possibly native 
    fid = fopen(filename,'r'); 
    head = fread(fid,4,'ushort'); % indeed mac native reread it 
    dim = head(4);                % figure out the dimensionality
    fseek(fid,0,'eof');
    f_size = (ftell(fid)-24)/4.0;
    stat = fseek(fid,0,-1);       % takes you back to the beginning of
    
    if dim == 2                   % the file
       head = fread(fid,12,'ushort');  % read the entire header
       n = head(6);
       m = f_size/n; 
       stat = fseek(fid,0,-1);            % back again to the front
       stat = fseek(fid,24,-1);
       a = fread(fid,[n,m],'single');    % finally read in the data
       a = a';
    elseif dim == 3
       head = char(fread(fid,12,'ushort'))';    % entire header
       n = head(6);                      % main dimensions
       m = head(8);
       l = head(10);                     % loop dimension
       % complete this after figuring out how to make 3d arrays
    end
    
  else       % not made on the mac

    head = fread(fid,20,'long');
    dim = head(2);
    if dim == 2
      n = head(3);
      m = head(4);
      stat = fseek(fid,0,-1);
      stat = fseek(fid,24,-1);
      a = fread(fid,[n,m],'single');
      a = a';
    elseif dim == 3
      n = head(6);                      % main dimensions
      m = head(8);
      l = head(10);     
      % complete this after figuring out how to make 3d arrays
    end
  end 

fclose(fid);



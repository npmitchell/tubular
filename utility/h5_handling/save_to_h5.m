function save_to_h5(h5fn, name, dataset, msg)
%SAVE_TO_H5 Wrapper for h5write to save safely.
%   Write the data in dataset to a dataset named name in h5fn. 
%   If a dataset already exists in h5fn, displays msg.

try
    h5create(h5fn, name, size(dataset)) ;
catch
    disp(msg)
end
try
    h5write(h5fn, name, dataset) ;
catch
    % The dataset was presumably the wrong size. Handle this case here
    fileattrib(h5fn,'+w');
    fid = H5F.open(h5fn, 'H5F_ACC_RDWR', 'H5P_DEFAULT');
    dset_id = H5D.open(fid, name);
    % NOTE: H5D.set_extent uses C-style indexing: [columns,rows]
    H5D.set_extent(dset_id, size(dataset));
    H5D.close(dset_id);
    H5F.close(fid);

    % Now try again
    h5write(h5fn, name, dataset) ;
    
    % todo: this doesn't work if contiguous storage...
    
end


function compile_mesh_self_intersection_3d
    mex -v -O mesh_self_intersection_3d.cpp -I/usr/include:/usr/local/include -L/usr/lib:/usr/local/lib -lCGAL -lgmp -lboost_thread
end


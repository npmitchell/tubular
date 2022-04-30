function compile_isotropic_remeshing
    mex -v -O isotropic_remeshing.cpp -I/usr/include:/usr/local/include -L/usr/lib:/usr/local/lib -lCGAL -lgmp -lboost_thread
end

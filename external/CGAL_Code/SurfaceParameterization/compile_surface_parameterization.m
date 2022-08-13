function compile_surface_parameterization
    mex -v -O surface_parameterization.cpp -I/usr/include:/usr/local/include -L/usr/lib:/usr/local/lib -lCGAL -lgmp -lboost_thread
end
function compile_surface_parameterization
    % As of 2025, we had this:
    % mex -v -O surface_parameterization.cpp -I/usr/include:/usr/local/include -L/usr/lib:/usr/local/lib -lCGAL -lgmp -lboost_thread
    % 2025-05-01 Now we removed -lCGAL since CGAL is header-only.
    % Boris also added CXXFLAGS="\$CXXFLAGS -std=c++17 -DCGAL_EIGEN3_ENABLED" 
    % because he downloaded eigen3 and wanted to link it explicitly
    mex -v -O CXXFLAGS="\$CXXFLAGS -std=c++17 -DCGAL_EIGEN3_ENABLED" -I/usr/include/eigen3 surface_parameterization.cpp -I/usr/include -I/usr/local/include -L/usr/lib -L/usr/local/lib -L/usr/share -lgmp -lboost_thread -lboost_system
end
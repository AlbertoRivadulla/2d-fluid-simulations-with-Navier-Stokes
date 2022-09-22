#ifndef NAVIER_STOKES_SOLVER_H
#define NAVIER_STOKES_SOLVER_H

#include <iostream>
#include <string>

#include "boundaryConditions.h"

class NavierStokesSolver
{
    public:
        // Constructor
        // Arguments: 
        //  - Dimensions
        //  - Boundary conditions
        //  - File to output the results
        //  - Seconds to run the simulation
        NavierStokesSolver( const int& Nx, const int& Ny, 
                            BoundaryConditions boundaryConditions,
                            const std::string& outFile,
                            const float& timeSim );

        // Run the simulation
        void solve();

    private:

};

#endif

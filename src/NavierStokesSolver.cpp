#include "NavierStokesSolver.h"

//=======================================================
//
// NavierStokesSolver class
//
//=======================================================

// Constructor
// Arguments: 
//  - Dimensions
//  - Boundary conditions
//  - File to output the results
//  - Seconds to run the simulation
NavierStokesSolver::NavierStokesSolver( const int& Nx, const int& Ny, 
                                        BoundaryConditions boundaryConditions,
                                        const std::string& outFile,
                                        const float& timeSim )
{
    // Check if the boundary conditions are correct
    if ( !boundaryConditions.check() )
    {
        std::cout << "Some boundary conditions for the walls are missing!\n";
        return;
    }

    // Initialize matrices

    //
    //
    //

}

// Run the simulation
void NavierStokesSolver::solve()
{
    // Iterate over time

        //
        //
        //

        // In each iteration, save the results to the file

}






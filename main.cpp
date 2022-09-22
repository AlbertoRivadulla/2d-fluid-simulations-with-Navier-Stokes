#include <iostream>
#include <string>
#include <vector>

#include "boundaryConditions.h"
#include "NavierStokesSolver.h"

int main()
{
    //---------------------------------------------------------------
    // Define the problem

    // Dimensions of the grid
    int Nx = 10;
    int Ny = 10;

    // Time to run the simulation
    double timeSim = 10;

    // Setup the boundary conditions
    BoundaryConditions boundaryConditions;
    boundaryConditions.addInflowWall( N, { 1, 100 } );
    boundaryConditions.addOutflowWall( S );
    boundaryConditions.addSolidWall( E );
    boundaryConditions.addSolidWall( W );

    //---------------------------------------------------------------
    // Run the simulation

    // Initialize the solver
    // Arguments: 
    //  - Dimensions
    //  - Boundary conditions
    //  - File to output the results
    //  - Seconds to run the simulation
    NavierStokesSolver solver( Nx, Ny, boundaryConditions, "out.txt", timeSim );

    // Run the solver

    //
    //
    //

    //---------------------------------------------------------------
    // Save the results to a video

    //
    //
    //

    return 0;
}

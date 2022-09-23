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
    double gridSize = 0.01;

    // Reynolds number
    double reynolds = 100.;

    // Time to run the simulation
    double timeSim = 2;
    // Parameter tau for the time step 
    double tau = 0.5;

    // Setup the boundary conditions
    BoundaryConditions boundaryConditions;
    boundaryConditions.addInflowWall( N, { 1, 100 } );
    boundaryConditions.addOutflowWall( S );
    boundaryConditions.addSolidWall( E );
    boundaryConditions.addSolidWall( W );

    // Add some solid objects in the domain 
    
    //

    //---------------------------------------------------------------
    // Run the simulation

    // Initialize the solver
    // Arguments: 
    //  - Dimensions
    //  - Boundary conditions
    //  - File to output the results
    //  - Seconds to run the simulation
    NavierStokesSolver solver( Nx, Ny, gridSize, reynolds, boundaryConditions, timeSim, tau );

    // Run the solver
    solver.solve( "out.txt" );

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

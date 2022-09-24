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
    int Nx = 5;
    int Ny = 5;
    double gridSize = 0.01;

    // Reynolds number
    double reynolds = 100.;

    // External force ( gravity )
    Vec2d force( 0., 0. );

    // Time to run the simulation
    double timeSim = 2;
    // double timeSim = 2./30.;
    // Parameter tau for the time step 
    double tau = 0.1;

    // Setup the boundary conditions
    BoundaryConditions boundaryConditions;
    // boundaryConditions.addInflowWall( N, { 1, 0 } );
    // boundaryConditions.addOutflowWall( S );
    // boundaryConditions.addSolidWall( E );
    // boundaryConditions.addSolidWall( W );
    boundaryConditions.addWallBC( N, INFLOW, { 1., 0. } );
    boundaryConditions.addWallBC( S, SOLID );
    boundaryConditions.addWallBC( E, SOLID );
    boundaryConditions.addWallBC( W, SOLID );

    // Add some solid objects in the domain 
    
    //
    //
    //

    //---------------------------------------------------------------
    // Run the simulation

    // Initialize the solver
    // Arguments: 
    //  - Dimensions
    //  - Boundary conditions
    //  - File to output the results
    //  - Seconds to run the simulation
    NavierStokesSolver solver( Nx, Ny, gridSize, reynolds, force, 
                               boundaryConditions, timeSim, tau );

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

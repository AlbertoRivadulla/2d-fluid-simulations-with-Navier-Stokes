#include <iostream>
#include <string>
#include <vector>

#include "boundaryConditions.h"
#include "NavierStokesSolver.h"
#include "videoOutput.h"

int main()
{
    //---------------------------------------------------------------
    // Define the problem

    // Dimensions of the grid
    int Nx = 200;
    int Ny = 150;
    double gridSize = 0.01;

    // Time to run the simulation
    double timeSim = 5.;

    // Parameter tau for the time step 
    double tau = 0.2;

    // Reynolds number
    double reynolds = 1000.;

    // External force ( gravity )
    Vec2d force( 0., 0. );

    // Setup the boundary conditions
    BoundaryConditions boundaryConditions;

    // boundaryConditions.addWallBC( N, INFLOW, { 0.0001, 0. } );
    // boundaryConditions.addWallBC( S, SOLID );
    // boundaryConditions.addWallBC( E, SOLID );
    // boundaryConditions.addWallBC( W, SOLID );

    boundaryConditions.addWallBC( W, INFLOW, { 0.000001, 0. } );
    boundaryConditions.addWallBC( N, SOLID );
    boundaryConditions.addWallBC( S, SOLID );
    // boundaryConditions.addWallBC( E, SOLID );

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

    //---------------------------------------------------------------
    // Save the results to a video

    // Create a video from the data written to a text file
    makeVideoFromDataFile( "out.txt", "video.mp4", 5 );

    return 0;
}

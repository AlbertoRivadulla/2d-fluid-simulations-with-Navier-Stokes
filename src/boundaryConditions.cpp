#include "boundaryConditions.h"


// Constructor
BoundaryConditions::BoundaryConditions()
{

}

// Add a solid wall
void BoundaryConditions::addSolidWall( Direction dir )
{
    mSolidWalls.push_back( dir );

    // Add the direction to the set of walls with BCs
    mWallsWithBCs.insert( dir );
}
// Add an inflow wall 
void BoundaryConditions::addInflowWall( Direction dir, const Vec2d& vel )
{
    mInflowWalls.push_back( std::pair<Direction, Vec2d>( dir, vel ) );

    // Add the direction to the set of walls with BCs
    mWallsWithBCs.insert( dir );
}
// Add an outflow wall 
void BoundaryConditions::addOutflowWall( Direction dir )
{
    mOutflowWalls.push_back( dir );

    // Add the direction to the set of walls with BCs
    mWallsWithBCs.insert( dir );
}

// Method to check if the boundary conditions are valid.
// They are valid if all the walls have some sort of BC.
bool BoundaryConditions::check()
{
    if ( mWallsWithBCs.size() < 4 )
        return false;
    return true;
}

// Method to apply the boundary conditions on an array of velocities
void BoundaryConditions::applyBCs( double* u, double* v )
{

    //
    //
    //

}


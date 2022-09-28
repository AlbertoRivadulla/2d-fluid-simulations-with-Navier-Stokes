#include "boundaryConditions.h"


// Constructor
BoundaryConditions::BoundaryConditions() :
    mNWall ( N, OUTFLOW, { 0., 0. } ),
    mSWall ( S, OUTFLOW, { 0., 0. } ),
    mEWall ( E, OUTFLOW, { 0., 0. } ),
    mWWall ( W, OUTFLOW, { 0., 0. } )
{
}

// Add boundary conditions to the walls
void BoundaryConditions::addWallBC( Direction dirVal, TypeBC typeVal, Vec2d velocityVal )
{
    switch ( dirVal )
    {
        case N:
            mNWall.type = typeVal;
            mNWall.velocity = velocityVal;
            break;
        case S:
            mSWall.type = typeVal;
            mSWall.velocity = velocityVal;
            break;
        case E:
            mEWall.type = typeVal;
            mEWall.velocity = velocityVal;
            break;
        case W:
            mWWall.type = typeVal;
            mWWall.velocity = velocityVal;
            break;
    }
}

// // add a solid wall
// void BoundaryConditions::addSolidWall( Direction dir )
// {
//     mSolidWalls.push_back( dir );
//
//     // Add the direction to the set of walls with BCs
//     mWallsWithBCs.insert( dir );
// }
// // Add an inflow wall 
// void BoundaryConditions::addInflowWall( Direction dir, const Vec2d& vel )
// {
//     mInflowWalls.push_back( std::pair<Direction, Vec2d>( dir, vel ) );
//
//     // Add the direction to the set of walls with BCs
//     mWallsWithBCs.insert( dir );
// }
// // Add an outflow wall 
// void BoundaryConditions::addOutflowWall( Direction dir )
// {
//     mOutflowWalls.push_back( dir );
//
//     // Add the direction to the set of walls with BCs
//     mWallsWithBCs.insert( dir );
// }

// // Method to check if the boundary conditions are valid.
// // They are valid if all the walls have some sort of BC.
// bool BoundaryConditions::check()
// {
//     if ( mWallsWithBCs.size() < 4 )
//         return false;
//     return true;
// }

// Method to apply the boundary conditions on an array of velocities
void BoundaryConditions::applyBCs( double* u, double* v, const int& Nx, const int& Ny )
{
    // ------------------------------------------------------------
    // Boundary conditions in the walls
    int Nyp2 = Ny + 2;

    // Top wall
    if ( mNWall.type != OUTFLOW )
    {
        for ( int i = 0; i < Nx + 2; ++i )
        {
            u[ i*Nyp2 + Ny+1 ] = 2. * mNWall.velocity.x - u[ i*Nyp2 + Ny ];
            v[ i*Nyp2 + Ny ]   = mNWall.velocity.y;
            v[ i*Nyp2 + Ny+1 ] = 2. * mNWall.velocity.y - v[ i*Nyp2 + Ny-1 ];
        }
    }
    else
    {
        // std::cout << "algoN\n";
        for ( int i = 0; i < Nx + 2; ++i )
        {
            u[ i*Nyp2 + Ny+1 ] = 2. * u[ i*Nyp2 + Ny ] - u[ i*Nyp2 + Ny-1 ];
            v[ i*Nyp2 + Ny ] = 2. * v[ i*Nyp2 + Ny-1 ] - v[ i*Nyp2 + Ny-2 ];
            v[ i*Nyp2 + Ny+1 ] = 2. * v[ i*Nyp2 + Ny ] - v[ i*Nyp2 + Ny-1 ];
        }
    }

    // Bottom wall
    if ( mSWall.type != OUTFLOW )
    {
        for ( int i = 0; i < Nx + 2; ++i )
        {
            u[ i*Nyp2 ] = 2. * mSWall.velocity.x - u[ i*Nyp2 + 1 ];
            v[ i*Nyp2 ] = mSWall.velocity.y; 
        }
    }
    else
    {
        // std::cout << "algoS\n";
        for ( int i = 0; i < Nx + 2; ++i )
        {
            u[ i*Nyp2 ] = 2. * u[ i*Nyp2 + 1 ] - u[ i*Nyp2 + 2 ];
            v[ i*Nyp2 ] = 2. * v[ i*Nyp2 + 1 ] - v[ i*Nyp2 + 2 ];
        }
    }

    // Left wall 
    if ( mWWall.type != OUTFLOW )
    {
        for ( int j = 0; j < Ny + 2; ++j )
        {
            u[ j ] = mWWall.velocity.x;
            v[ j ] = 2. * mWWall.velocity.y - v[ Nyp2 + j ];
        }
    }
    else
    {
        for ( int j = 0; j < Ny + 2; ++j )
        {
            u[ j ] = 2. * u[ Nyp2 + j ] - u[ 2*Nyp2 + j ];
            v[ j ] = 2. * v[ Nyp2 + j ] - v[ 2*Nyp2 + j ];
        }
    }

    // Right wall
    if ( mEWall.type != OUTFLOW )
    {
        for ( int j = 0; j < Ny + 2; ++j )
        {
            u[ Nx*Nyp2 + j ]     = mEWall.velocity.x;
            u[ (Nx+1)*Nyp2 + j ] = 2. * mEWall.velocity.x - u[ (Nx-1)*Nyp2 + j ];
            v[ (Nx+1)*Nyp2 + j ] = 2. * mEWall.velocity.y - v[ Nx*(Ny+1) + j ];
        }
    }
    else
    {
        for ( int j = 0; j < Ny + 2; ++j )
        {
            u[ Nx*Nyp2 + j ] = 2. * u[ (Nx-1)*Nyp2 + j ] - u[ (Nx-2)*Nyp2 + j ];
            u[ (Nx+1)*Nyp2 + j ] = 2. * u[ Nx*Nyp2 + j ] - u[ (Nx-1)*Nyp2 + j ];
            v[ (Nx+1)*Nyp2 + j ] = 2. * v[ Nx*Nyp2 + j ] - v[ (Nx-1)*Nyp2 + j ];
        }
    }


    // ------------------------------------------------------------
    // Boundary conditions due to objects inside the domain

    //
    //
    //

}


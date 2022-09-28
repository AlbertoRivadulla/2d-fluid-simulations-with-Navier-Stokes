#include "boundaryConditions.h"


// Constructor
BoundaryConditions::BoundaryConditions( const int& Nx, const int& Ny ) :
    mNx { Nx }, mNy { Ny },
    mNWall ( N, OUTFLOW, { 0., 0. } ),
    mSWall ( S, OUTFLOW, { 0., 0. } ),
    mEWall ( E, OUTFLOW, { 0., 0. } ),
    mWWall ( W, OUTFLOW, { 0., 0. } ),
    // Initialize the domain map with ones, which correspond to fluid cells
    mDomainMap ( Nx * Ny, 1 )
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

void BoundaryConditions::addObstacleCircle( const int& x, const int& y, const double& radius )
{
    // Draw a circle on the domain wall 
    for ( int i = 0; i < mNx; ++i )
        for ( int j = 0; j < mNy; ++j )
            if ( (i-x)*(i-x) + (j-y)*(j-y) <= radius*radius )
                mDomainMap[ i*mNy + j ] = 0;

    // Add the boundary to the list of boundary conditions

    //
    //
    //
}
// void addObstacleRectangle( const int& x, const int& y, const double& width, const double& height )
// {
//
// }


// Save the map of the domain to a file stream
void BoundaryConditions::drawDomainMap( std::ofstream& fileStream )
{
    for ( int i = 0; i < mNx*mNy; ++i )
    {
        if ( mDomainMap[ i ] == 0 )
            fileStream << '0';
        else
            fileStream << '1';
    }
    fileStream << '\n';
}

// Method to check if a pair of indices correspond to a fluid cell
bool BoundaryConditions::isFluidCell( const int& i, const int& j )
{
    return mDomainMap[ i*mNy + j ];
}

// Method to apply the boundary conditions on an array of velocities
void BoundaryConditions::applyBCs( double* u, double* v )
{
    // Boundary conditions in the walls
    applyWallsBCs( u, v );

    // Boundary conditions due to objects inside the domain
    applyObstaclesBCs( u, v );
}

// Apply boundary conditions due to the walls
void BoundaryConditions::applyWallsBCs( double* u, double* v )
{
    int Nyp2 = mNy + 2;

    // Top wall
    if ( mNWall.type != OUTFLOW )
    {
        for ( int i = 0; i < mNx + 2; ++i )
        {
            u[ i*Nyp2 + mNy+1 ] = 2. * mNWall.velocity.x - u[ i*Nyp2 + mNy ];
            v[ i*Nyp2 + mNy ]   = mNWall.velocity.y;
            v[ i*Nyp2 + mNy+1 ] = 2. * mNWall.velocity.y - v[ i*Nyp2 + mNy-1 ];
        }
    }
    else
    {
        for ( int i = 0; i < mNx + 2; ++i )
        {
            u[ i*Nyp2 + mNy+1 ] = 2. * u[ i*Nyp2 + mNy ] - u[ i*Nyp2 + mNy-1 ];
            v[ i*Nyp2 + mNy ] = 2. * v[ i*Nyp2 + mNy-1 ] - v[ i*Nyp2 + mNy-2 ];
            v[ i*Nyp2 + mNy+1 ] = 2. * v[ i*Nyp2 + mNy ] - v[ i*Nyp2 + mNy-1 ];
        }
    }

    // Bottom wall
    if ( mSWall.type != OUTFLOW )
    {
        for ( int i = 0; i < mNx + 2; ++i )
        {
            u[ i*Nyp2 ] = 2. * mSWall.velocity.x - u[ i*Nyp2 + 1 ];
            v[ i*Nyp2 ] = mSWall.velocity.y; 
        }
    }
    else
    {
        for ( int i = 0; i < mNx + 2; ++i )
        {
            u[ i*Nyp2 ] = 2. * u[ i*Nyp2 + 1 ] - u[ i*Nyp2 + 2 ];
            v[ i*Nyp2 ] = 2. * v[ i*Nyp2 + 1 ] - v[ i*Nyp2 + 2 ];
        }
    }

    // Left wall 
    if ( mWWall.type != OUTFLOW )
    {
        for ( int j = 0; j < mNy + 2; ++j )
        {
            u[ j ] = mWWall.velocity.x;
            v[ j ] = 2. * mWWall.velocity.y - v[ Nyp2 + j ];
        }
    }
    else
    {
        for ( int j = 0; j < mNy + 2; ++j )
        {
            u[ j ] = 2. * u[ Nyp2 + j ] - u[ 2*Nyp2 + j ];
            v[ j ] = 2. * v[ Nyp2 + j ] - v[ 2*Nyp2 + j ];
        }
    }

    // Right wall
    if ( mEWall.type != OUTFLOW )
    {
        for ( int j = 0; j < mNy + 2; ++j )
        {
            u[ mNx*Nyp2 + j ]     = mEWall.velocity.x;
            u[ (mNx+1)*Nyp2 + j ] = 2. * mEWall.velocity.x - u[ (mNx-1)*Nyp2 + j ];
            v[ (mNx+1)*Nyp2 + j ] = 2. * mEWall.velocity.y - v[ mNx*(mNy+1) + j ];
        }
    }
    else
    {
        for ( int j = 0; j < mNy + 2; ++j )
        {
            u[ mNx*Nyp2 + j ] = 2. * u[ (mNx-1)*Nyp2 + j ] - u[ (mNx-2)*Nyp2 + j ];
            u[ (mNx+1)*Nyp2 + j ] = 2. * u[ mNx*Nyp2 + j ] - u[ (mNx-1)*Nyp2 + j ];
            v[ (mNx+1)*Nyp2 + j ] = 2. * v[ mNx*Nyp2 + j ] - v[ (mNx-1)*Nyp2 + j ];
        }
    }
}

// Apply boundary conditions due to obstacles
void BoundaryConditions::applyObstaclesBCs( double* u, double* v )
{

    //
    //
    //

}

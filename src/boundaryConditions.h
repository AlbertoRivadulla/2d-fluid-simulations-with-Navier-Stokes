#ifndef BOUNDARY_CONDITIONS_H
#define BOUNDARY_CONDITIONS_H

#include <iostream>
#include <vector>
#include <set>
#include <utility> // For std::pair
#include <fstream>

// 2-dimensional vector for the velocities
struct Vec2d
{
    double x;
    double y;

    // Constructor
    Vec2d() : x { 0. }, y { 0. }
    {}
    Vec2d ( const double& xVal, const double& yVal )
    {
        x = xVal;
        y = yVal;
    }
};

// Four directions as cardinal points
enum Direction
{
    N,
    S,
    E,
    W,
};

// Types of boundary conditions
enum TypeBC
{
    INFLOW,
    OUTFLOW,
    SOLID,
};

// Boundary condition in a wall 
struct WallBC
{
    // Attributes
    Direction dir;
    Vec2d velocity;
    TypeBC type;

    // Constructor
    WallBC ( Direction dirVal, TypeBC typeVal, Vec2d velocityVal ) :
        dir { dirVal }, velocity { velocityVal }, type { typeVal }
    {
    }
};

// Struct for the boundary conditions
class BoundaryConditions
{
    public:
        // Constructor
        // By default the 4 wall will be OUTFLOW walls
        BoundaryConditions( const int& Nx, const int& Ny );

        // Add boundary conditions to the walls
        void addWallBC( Direction dirVal, TypeBC typeVal, Vec2d velocityVal = { 0., 0. } );

        // Methods to add different obstacles
        void addObstacleCircle( const int& x, const int& y, const double& radius );
        // void addObstacleRectangle( const int& x, const int& y, const double& width, const double& height );

        // Save the map of the domain to a file stream
        void drawDomainMap( std::ofstream& fileStream );

        // Method to check if a pair of indices correspond to a fluid cell
        bool isFluidCell( const int& i, const int& j );

        // Method to apply the boundary conditions on an array of velocities
        void applyBCs( double* u, double* v );

    private:
        // Dimensions of the domain
        int mNx;
        int mNy;

        // Map of the domain
        //      0: obstacle cell
        //      1: fluid cell
        // Walls are not considered here, but outside
        // unsigned int* mDomainMap;
        std::vector<unsigned int> mDomainMap;

        // Boundary conditions for the 4 walls
        WallBC mNWall;
        WallBC mSWall;
        WallBC mEWall;
        WallBC mWWall;

        // Apply boundary conditions due to the walls
        void applyWallsBCs( double* u, double* v );

        // Apply boundary conditions due to obstacles
        void applyObstaclesBCs( double* u, double* v );
};

#endif

#ifndef BOUNDARY_CONDITIONS_H
#define BOUNDARY_CONDITIONS_H

#include <iostream>
#include <vector>
#include <set>
#include <utility> // For std::pair

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
        BoundaryConditions();

        // Add boundary conditions to the walls
        void addWallBC( Direction dirVal, TypeBC typeVal, Vec2d velocityVal = { 0., 0. } );

        // // Add a solid wall
        // void addSolidWall( Direction dir );
        // // Add an inflow wall 
        // void addInflowWall( Direction dir, const Vec2d& vel );
        // // Add an outflow wall 
        // void addOutflowWall( Direction dir );

        // // Method to check if the boundary conditions are valid.
        // // They are valid if all the walls have some sort of BC.
        // bool check();

        // Method to apply the boundary conditions on an array of velocities
        void applyBCs( double* u, double* v, const int& Nx, const int& Ny );

    private:
        // // List of walls with boundary conditions
        // std::set<Direction> mWallsWithBCs;
        //
        // // Solid walls
        // std::vector<Direction> mSolidWalls;
        //
        // // Inflow walls
        // std::vector<std::pair<Direction, Vec2d>> mInflowWalls;
        //
        // // Outflow walls 
        // std::vector<Direction> mOutflowWalls;

        // Boundary conditions for the 4 walls
        WallBC mNWall;
        WallBC mSWall;
        WallBC mEWall;
        WallBC mWWall;
};

#endif

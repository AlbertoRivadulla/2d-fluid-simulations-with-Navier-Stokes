#ifndef BOUNDARY_CONDITIONS_H
#define BOUNDARY_CONDITIONS_H

#include <iostream>
#include <vector>
#include <set>
#include <utility> // For std::pair

// Four directions as cardinal points
enum Direction
{
    N,
    S,
    E,
    W,
};

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

// Struct for the boundary conditions
class BoundaryConditions
{
    public:
        // Constructor
        BoundaryConditions();

        // Add a solid wall
        void addSolidWall( Direction dir );
        // Add an inflow wall 
        void addInflowWall( Direction dir, const Vec2d& vel );
        // Add an outflow wall 
        void addOutflowWall( Direction dir );

        // Method to check if the boundary conditions are valid.
        // They are valid if all the walls have some sort of BC.
        bool check();

        // Method to apply the boundary conditions on an array of velocities
        void applyBCs( double* u, double* v );

    private:
        // List of walls with boundary conditions
        std::set<Direction> mWallsWithBCs;

        // Solid walls
        std::vector<Direction> mSolidWalls;

        // Inflow walls
        std::vector<std::pair<Direction, Vec2d>> mInflowWalls;

        // Outflow walls 
        std::vector<Direction> mOutflowWalls;
};



#endif

#ifndef NAVIER_STOKES_SOLVER_H
#define NAVIER_STOKES_SOLVER_H

#include <iostream>
#include <string>
#include <limits>
#include <algorithm>

#include "boundaryConditions.h"

class NavierStokesSolver
{
    public:
        // Constructor
        // Arguments: 
        //  - Dimensions
        //  - Boundary conditions
        //  - File to output the results
        //  - Seconds to run the simulation
        NavierStokesSolver( const int& Nx, const int& Ny, const double& gridSize,
                            const double& reynolds,
                            BoundaryConditions boundaryConditions,
                            const double& timeSim,
                            const double& tau = 0.5 );

        // Methods to set other paramters
        void setFrameRate( double fps );
        void setOmegaSOR( double omega );

        // Set initial values
        void setInitialValues( double u0, double v0, double p0 );

        // Run the simulation
        void solve( std::string outFile );

    private:

        // --------------------------------------------------------------
        // Paramters of the simulation

        // Size of the grid
        int mNx, mNy;
        double mH;
        // Reynolds number
        double mReynolds;
        // Pointer to the boundary conditions
        BoundaryConditions* mBoundaryConditions;
        // Time to run the simulation
        double mTimeSim;
        // Parameter tau for the time step 
        double mTau;
        // Frame time 
        double mFrameTime;
        // Parameter omega for SOR (successive over-relaxation)
        double mOmega;

        // --------------------------------------------------------------
        // Parameters used during the simulation
        double mDeltaT;
        double mThisFrameTime;
        double mTotalTime;

        int mTotalFrames;
        int mFrameCount;

        // --------------------------------------------------------------
        // Matrices

        // Velocities in the x (u) and y (v) directions
        double* mU;
        double* mV;

        // Matrix R in the x and y directions
        // Two copies, to hold the new and the old ones
        double* mRx1;
        double* mRy1;
        double* mRx2;
        double* mRy2;
        // Pointers to the new and old matrices
        double* mRx;
        double* mRy;
        double* mRxOld;
        double* mRyOld;

        // Tentative velocities u^p and v^p
        double* mUp;
        double* mVp;

        // Auxiliary rhs matrix
        double* mRHS;

        // Pressure (actually, this is the pseudo-pressure, equal to p*deltat)
        // Two copies, to hold the new and the old ones
        double* mP1;
        double* mP2;
        // Pointers to the new and old matrices
        double* mP;
        double* mPOld;

        // --------------------------------------------------------------
        // Private methods

        // Method to initialize the matrices
        void initializeMatrices();

        // Compute the delta of time
        void computeDeltaTime();

        // Compute the tentative velocity
        void computeTentativeVelocity();

        // Compute the RHS for the Poisson eq. for the pressure 
        void computeRHSPoisson();

        // Solve the Poisson eq. for the pressure 
        void solvePoisson();

        // Update the velocity
        void updateVelocity();

        // Save the results to the file
        void saveStepToFile( const std::string& outFile );
};

#endif

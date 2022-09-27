#ifndef NAVIER_STOKES_SOLVER_H
#define NAVIER_STOKES_SOLVER_H

#include <iostream>
#include <string>
#include <limits>
#include <algorithm>
#include <cmath>
#include <chrono>
#include <fstream>

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
                            const Vec2d& force,
                            BoundaryConditions boundaryConditions,
                            const double& timeSim,
                            const double& tau = 0.5 );

        // Destructor
        ~NavierStokesSolver();

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
        double mH, mHInv, mHSq, mHSqInv;
        // Size of the grid with the extra rows and columns
        int mNxp2, mNyp2;
        // Reynolds number
        double mReynolds, mReynoldsInv;
        // Gravity and external forces
        double mGx, mGy;
        // Pointer to the boundary conditions
        BoundaryConditions* mBoundaryConditions;
        // Time to run the simulation
        double mTimeSim;
        // Parameter tau for the time step 
        double mTau;
        // Frame time 
        double mFrameTime;
        // Parameter omega for SOR (successive over-relaxation),
        // threshold for the L2-norm of the residual squared and
        // maximum number of iterations
        double mOmega;
        double mThreshold;
        double mIterMax;

        // --------------------------------------------------------------
        // Parameters used during the simulation
        double mDeltaT;
        double mThisFrameTime;
        double mTotalTime;

        int mTotalFrames;
        int mFrameCount;

        // Boolean variable to know if the simulation is good
        bool mGoodSimulation;

        // Scale of the velocity, used to scale it in the video
        double mVelocityScale;

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

        // Pressure
        double* mP;

        // --------------------------------------------------------------
        // Variables used for the output

        // File to write
        std::ofstream mOutputFile;

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

        // Initialize the file to save the results
        void initOutputFile( const std::string& outFileName );

        // Save the results to the file
        void saveStepToFile();

        // Write the velocity scale at the beginning of the file 
        void writeVelocityScale();
};

#endif

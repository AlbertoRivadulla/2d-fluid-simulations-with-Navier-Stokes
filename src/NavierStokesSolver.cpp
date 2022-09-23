#include "NavierStokesSolver.h"

//=======================================================
//
// NavierStokesSolver class
//
//=======================================================

// ------------------------------------------------------
// Initializer methods

// Constructor
// Arguments: 
//  - Dimensions
//  - Boundary conditions
//  - File to output the results
//  - Seconds to run the simulation
NavierStokesSolver::NavierStokesSolver( const int& Nx, const int& Ny, 
                                        const double& gridSize,
                                        const double& reynolds,
                                        BoundaryConditions boundaryConditions,
                                        const double& timeSim,
                                        const double& tau )
{
    // Check if the boundary conditions are correct
    if ( !boundaryConditions.check() )
    {
        std::cout << "Some boundary conditions for the walls are missing!\n";
        return;
    }

    // Store the given parameters
    mNx = Nx;
    mNy = Ny;
    mH = gridSize;
    mReynolds = reynolds;
    mBoundaryConditions = &boundaryConditions;
    mTimeSim = timeSim;
    mTau = tau;

    // Initialize other parameters
    // Frame time - default 30 FPS
    mFrameTime = 1. / 30.;

    // Omega - parameter for SOR (Successve Over-Relaxation)
    mOmega = 1.7;

    // Parameters of the time of the simulation
    mDeltaT = 0.;
    mThisFrameTime = 0.;
    mTotalTime = 0.;
    mTotalFrames = (int)(mTimeSim / mFrameTime);
    mFrameCount = 0;

    // Initialize matrices
    initializeMatrices();
    setInitialValues( 0., 0., 0. );
}

// Methods to set other paramters
void NavierStokesSolver::setFrameRate( double fps )
{
    mFrameTime = 1. / fps;
    mTotalFrames = (int)( fps * mTotalTime );
}
void NavierStokesSolver::setOmegaSOR( double omega )
{
    mOmega = omega;
}

// Set initial values
void NavierStokesSolver::setInitialValues( double u0, double v0, double p0 )
{
    for ( int i = 0; i < mNx + 2; ++i )
    {
        for ( int j = 0; j < mNy + 2; ++j )
        {
            mU[ i*mNx + j ] = u0;
            mV[ i*mNx + j ] = v0;
            mP[ i*mNx + j ] = p0;
        }
    }
}

// ------------------------------------------------------
// Main loop

// Run the simulation
void NavierStokesSolver::solve( std::string outFile )
{
    // Check if the boundary conditions are correct
    if ( !mBoundaryConditions->check() )
    {
        std::cout << "Can't run the simulation with some boundary conditions missing!\n";
        return;
    }
    // Iterate over time
    while ( mTotalTime < mTimeSim )
    {
        // Compute the delta of time of the frame
        computeDeltaTime();

        mThisFrameTime += mDeltaT;
        mTotalTime += mDeltaT;

        // Compute the tentative velocity
        computeTentativeVelocity();

        // Compute the RHS for the Poisson eq. for the pressure 
        computeRHSPoisson();

        // Solve the Poisson eq. for the pressure 
        solvePoisson();

        // Update the velocity
        updateVelocity();

        // Apply the boundary conditions
        mBoundaryConditions->applyBCs( mU, mV );

        // If the time elapsed is more than a frame, save the current state
        if ( mThisFrameTime >= mFrameTime )
        {
            mFrameCount += 1;
            mThisFrameTime -= mFrameTime;

            // Save the results to the file
            saveStepToFile( outFile );

            std::cout << "Frame " << mFrameCount << " of " << mTotalFrames << "," 
                      << mTotalTime << " seconds.\n";
        }

        // Swap old and new matrices for R and p
        std::swap( mRx, mRxOld );
        std::swap( mRy, mRyOld );
        std::swap( mP, mPOld );
        // double* tempRx;
        // double* tempRy;
        // tempRx = mRx;
        // tempRy = mRy;
        // mRx = mRxOld;
        // mRy = mRyOld;
        // mRxOld = tempRx;
        // mRyOld = tempRy;
        //
        // double* tempP;
        // tempP = mP;
        // mP = mPOld;
        // mPOld = tempP;
    }
}

// ------------------------------------------------------
// Private methods

// Method to initialize the matrices
void NavierStokesSolver::initializeMatrices()
{
    // Velocities in the x (u) and y (v) directions
    mU = new double [ (mNx + 2) * ( mNy + 2 ) ];
    mV = new double [ (mNx + 2) * ( mNy + 2 ) ];

    // Matrix R in the x and y directions
    // Two copies, to hold the new and the old ones
    mRx1 = new double [ mNx * mNy ];
    mRy1 = new double [ mNx * mNy ];
    mRx2 = new double [ mNx * mNy ];
    mRy2 = new double [ mNx * mNy ];
    // Pointers to the new and old matrices
    mRx = mRx1;
    mRy = mRy1;
    mRxOld = mRx2;
    mRyOld = mRy2;

    // Tentative velocities u^p and v^p
    mUp = new double [ (mNx + 2) * ( mNy + 2 ) ];
    mVp = new double [ (mNx + 2) * ( mNy + 2 ) ];

    // Auxiliary rhs matrix
    mRHS = new double [ mNx * mNy ];

    // Pressure (actually, this is the pseudo-pressure, equal to p*deltat)
    // Two copies, to hold the new and the old ones
    mP1 = new double [ (mNx + 2) * ( mNy + 2 ) ];
    mP2 = new double [ (mNx + 2) * ( mNy + 2 ) ];
    // Pointers to the new and old matrices
    mP = mP1;
    mPOld = mP2;
}

// Compute the delta of time
void NavierStokesSolver::computeDeltaTime()
{
    // Compute the maximum velocity in each direction
    double uMax = std::numeric_limits<double>::max();
    double vMax = std::numeric_limits<double>::max();
    for ( int i = 0; i < mNx + 2; ++i )
    {
        for ( int j = 0; j < mNy + 2; ++j )
        {
            if ( mU[ i*mNx + j ] > uMax )
                uMax = mU[ i*mNx + j ];
            if ( mV[ i*mNx + j ] > uMax )
                vMax = mV[ i*mNx + j ];
        }
    }
    // If the maximum velocity is zero, add a small offset to it 
    if ( uMax == 0. ) uMax = 0.000001;
    if ( vMax == 0. ) vMax = 0.000001;

    // Take the minimum of the three values 
    mDeltaT = std::min( mReynolds / 4. * mH, std::min( 1. / uMax, 1. / vMax ) );

    // Multiply the time step by tau, and by the step h that I left out in the previous line
    mDeltaT *= mTau * mH;
    
    // Check that the time step is smaller than the length of the frame
    if ( mDeltaT < mFrameTime )
        mDeltaT = mFrameTime;
}

// Compute the tentative velocity
void NavierStokesSolver::computeTentativeVelocity()
{
    
}

// Compute the RHS for the Poisson eq. for the pressure 
void NavierStokesSolver::computeRHSPoisson()
{
    
}

// Solve the Poisson eq. for the pressure 
void NavierStokesSolver::solvePoisson()
{
    
}

// Update the velocity
void NavierStokesSolver::updateVelocity()
{
    
}

// Save the results to the file
void NavierStokesSolver::saveStepToFile( const std::string& outFile )
{
    
}


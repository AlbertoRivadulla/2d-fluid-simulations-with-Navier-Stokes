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
                                        const Vec2d& force,
                                        BoundaryConditions boundaryConditions,
                                        const double& timeSim,
                                        const double& tau )
{
    // // Check if the boundary conditions are correct
    // if ( !boundaryConditions.check() )
    // {
    //     std::cout << "Some boundary conditions for the walls are missing!\n";
    //     return;
    // }

    // Store the given parameters
    mNx = Nx;
    mNy = Ny;
    mNxp2 = Nx + 2;
    mNyp2 = Ny + 2;
    mH = gridSize;
    mHInv = 1. / gridSize;
    mHSq = gridSize * gridSize;
    mHSqInv = 1. / mHSq;
    mReynolds = reynolds;
    mReynoldsInv = 1. / reynolds;
    mBoundaryConditions = &boundaryConditions;
    mTimeSim = timeSim;
    mTau = tau;
    mGx = force.x;
    mGy = force.y;

    // Initialize other parameters
    // Frame time - default 30 FPS
    mFrameTime = 1. / 30.;

    // Omega - parameter for SOR (Successve Over-Relaxation)
    mOmega = 1.7;
    // Threshold for the L2-norm of the residual squared
    mThreshold = 0.0001 * mNx * mNy;
    // Maximum number of iterations for the algorithm
    mIterMax = 100;

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
    for ( int i = 0; i < mNxp2; ++i )
    {
        for ( int j = 0; j < mNyp2; ++j )
        {
            mU[ i*mNyp2 + j ] = u0;
            mV[ i*mNyp2 + j ] = v0;
            mP[ i*mNyp2 + j ] = p0;
        }
    }
}

// ------------------------------------------------------
// Main loop

// Run the simulation
void NavierStokesSolver::solve( std::string outFile )
{
    // // Check if the boundary conditions are correct
    // if ( !mBoundaryConditions->check() )
    // {
    //     std::cout << "Can't run the simulation with some boundary conditions missing!\n";
    //     return;
    // }

    // Initialize the output file 
    initOutputFile( outFile );

    // Initialize the time counter
	auto t0 {std::chrono::high_resolution_clock::now()};

    // Apply the boundary conditions for the first time
    mBoundaryConditions->applyBCs( mU, mV, mNx, mNy );

    int count = 0;
    // Iterate over time
    while ( mTotalTime < mTimeSim )
    {
        // if ( count++ > 4 )
        //     break;
        // Compute the delta of time of the frame
        computeDeltaTime();

        // std::cout << "delta t: " << mDeltaT << std::endl;

        // Compute the tentative velocity
        computeTentativeVelocity();

        // Compute the RHS for the Poisson eq. for the pressure 
        computeRHSPoisson();

        // Solve the Poisson eq. for the pressure 
        solvePoisson();

        // Update the velocity
        updateVelocity();

        // Apply the boundary conditions
        mBoundaryConditions->applyBCs( mU, mV, mNx, mNy );

        // If the time elapsed is more than a frame, save the current state
        mThisFrameTime += mDeltaT;
        mTotalTime += mDeltaT;
        if ( mThisFrameTime >= mFrameTime )
        {
            mFrameCount += 1;
            mThisFrameTime -= mFrameTime;

            // Save the results to the file
            saveStepToFile( outFile );

            std::cout << "Frame " << mFrameCount << " of " << mTotalFrames << " --- " 
                      << mTotalTime << " seconds.\n";
        }

        // // Swap old and new matrices for R and p
        // std::swap( mRx, mRxOld );
        // std::swap( mRy, mRyOld );
        // std::swap( mP, mPOld );
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

    // Compute the total time of the simulation
	auto t1 {std::chrono::high_resolution_clock::now()};
	std::chrono::duration<double> duration {t1 - t0};
    std::cout << "\nSimulation done in " << duration.count() << " seconds\n";

    // Close the output file 
    closeOutputFile( outFile );
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
    mP = new double [ (mNx + 2) * (mNy + 2) ];
    // Two copies, to hold the new and the old ones
    // mP1 = new double [ (mNx + 2) * ( mNy + 2 ) ];
    // mP2 = new double [ (mNx + 2) * ( mNy + 2 ) ];
    // // Pointers to the new and old matrices
    // mP = mP1;
    // mPOld = mP2;
}

// Compute the delta of time
void NavierStokesSolver::computeDeltaTime()
{
    // Compute the maximum velocity in each direction
    double uMax = 0.;
    double vMax = 0.;
    for ( int i = 0; i < mNxp2; ++i )
    {
        for ( int j = 0; j < mNyp2; ++j )
        {
            if ( mU[ i*mNyp2 + j ] > uMax )
                uMax = mU[ i*mNyp2 + j ];
            if ( mV[ i*mNyp2 + j ] > uMax )
                vMax = mV[ i*mNyp2 + j ];
        }
    }
    // If the maximum velocity is zero, add a small offset to it 
    uMax = std::max( std::fabs( uMax ), 0.000001 );
    vMax = std::max( std::fabs( vMax ), 0.000001 );

    // Take the minimum of the three values 
    mDeltaT = std::min( mReynolds / 4. * mH, std::min( 1. / uMax, 1. / vMax ) );

    // Multiply the time step by tau, and by the step h that I left out in the previous line
    mDeltaT *= mTau * mH;
    
    // Check that the time step is smaller than the length of the frame
    mDeltaT = std::min( mDeltaT, mFrameTime );
}

// Compute the tentative velocity
void NavierStokesSolver::computeTentativeVelocity()
{
    // Compute Rx and Ry, and Up and Vp in the same loop
    for ( int i = 1; i < mNx + 1; ++i )
    {
        for ( int j = 1; j < mNy + 1; ++j )
        {
            // Compute indices
            unsigned int ij     = i * mNyp2 + j;
            unsigned int im1j   = ( i - 1 ) * mNyp2 + j;
            unsigned int ip1j   = ( i + 1 ) * mNyp2 + j;
            unsigned int ijm1   = i * mNyp2 + j - 1;
            unsigned int ijp1   = i * mNyp2 + j + 1;
            unsigned int ijForR = (i-1) * mNy + j-1;

            // Compute Rx and Ry
            mRx[ ijForR ] = - mU[ ij ] * ( mU[ ij ] - mU[ im1j ] ) * mHInv 
                            - mV[ ij ] * ( mU[ ij ] - mU[ ijm1 ] ) * mHInv 
                            + mReynoldsInv * mHSqInv * ( mU[ ip1j ] - 2.*mU[ ij ] + mU[im1j] + 
                                                         mU[ ijp1 ] - 2.*mU[ ij ] + mU[ijm1] )
                            + mGx;
            mRy[ ijForR ] = - mU[ ij ] * ( mV[ ij ] - mV[ im1j ] ) * mHInv 
                            - mV[ ij ] * ( mV[ ij ] - mV[ ijm1 ] ) * mHInv 
                            + mReynoldsInv * mHSqInv * ( mV[ ip1j ] - 2.*mV[ ij ] + mV[im1j] + 
                                                         mV[ ijp1 ] - 2.*mV[ ij ] + mV[ijm1] )
                            + mGy;

            // Compute Up and Vp 
            if ( mTotalTime == 0. )
            {
                mUp[ ij ] = mU[ ij ] + mDeltaT * mRx[ ijForR ];
                mVp[ ij ] = mV[ ij ] + mDeltaT * mRy[ ijForR ];
            }
            else
            {
                mUp[ ij ] = mU[ ij ] + mDeltaT * ( 1.5 * mRx[ ijForR ] - 0.5 * mRxOld[ ijForR ] );
                mVp[ ij ] = mV[ ij ] + mDeltaT * ( 1.5 * mRy[ ijForR ] - 0.5 * mRyOld[ ijForR ] );
            }

            // std::cout << mRx[ ijForR ] << ' ' << mRy[ ijForR ] << " /// ";
            // std::cout << mUp[ ij ] << ' ' << mVp[ ij ] << " === ";
        }
    }

    // // Compute Up and Vp
    // for ( int i = 1; i < mNx + 1; ++i )
    // {
    //     for ( int j = 1; j < mNy + 1; ++j )
    //     {
    //         unsigned int ij = i * mNyp2 + j;
    //         unsigned int ijForR = i*mNy + j;
    //         mUp[ ij ] = mU[ ij ] + mDeltaT * ( 1.5 * mRx[ ijForR ] - 0.5 * mRxOld[ ijForR ] );
    //         mVp[ ij ] = mV[ ij ] + mDeltaT * ( 1.5 * mRy[ ijForR ] - 0.5 * mRyOld[ ijForR ] );
    //     }
    // }

    // Set the values of the components up[0, j] and vp[i, 0] so that the derivatives
    // next to these boundaries are computed correctly when evaluating the RHS of
    // the Poisson equation 
    for ( int j = 1; j < mNy + 1; ++j )
        mUp[ j ] = 2. * mUp[ mNyp2 + j ] - mUp[ 2*mNyp2 + j ];
    for ( int i = 1; i < mNx + 1; ++i )
        mVp[ i*mNyp2 ] = 2. * mVp[ i*mNyp2 + 1 ] - mVp[ i*mNyp2 + 2 ];

    // Swap old and new matrices for Rx and Ry
    std::swap( mRx, mRxOld );
    std::swap( mRy, mRyOld );
}

// Compute the RHS for the Poisson eq. for the pressure 
void NavierStokesSolver::computeRHSPoisson()
{
    for ( int i = 1; i < mNx + 1; ++i )
    {
        for ( int j = 1; j < mNy + 1; ++j )
        {
            unsigned int ij       = i * mNyp2 + j;
            unsigned int im1j     = ( i - 1 ) * mNyp2 + j;
            unsigned int ijm1     = i * mNyp2 + j - 1;
            unsigned int ijForRHS = (i-1) * mNy + (j-1);
            // std::cout << mUp[ ij ] << ' ' << mUp[ im1j ] << ' ' << mVp[ ij ] << ' ' << mVp [ ij ] << " ----- ";
            mRHS[ ijForRHS ] = mHInv * ( mUp[ ij ] - mUp[ im1j ] + mVp[ ij ] - mVp[ ijm1 ] );

            // std::cout<< std::endl << "RHS: ";
            // std::cout << mRHS[ ijForRHS ] << " --- ";
        }
    }
}

// Solve the Poisson eq. for the pressure 
void NavierStokesSolver::solvePoisson()
{
    // Initialize the L2-norm of the residual
    double residualSq;
    double residualSqOld = 0.;

    // Set the old pressure to be the average of the squared pressure
    double pSqAcc = 0.;
    for ( int i = 1; i < mNx; ++i )
        for ( int j = 1; j < mNy; ++j )
            // pSqAcc += mPOld[ i*mNyp2 + j ] * mPOld[ i*mNyp2 + j ];
            pSqAcc += mP[ i*mNyp2 + j ] * mP[ i*mNyp2 + j ];
    pSqAcc /= mNx * mNy;
    pSqAcc = std::sqrt( pSqAcc );

    for ( int i = 1; i < mNx + 1; ++i )
        for ( int j = 1; j < mNy + 1; ++j )
            mP[ i*mNyp2 + j ] = pSqAcc;

    // std::cout << std::endl << "Initial pressure: ";
    // for ( int i = 1; i < mNx + 1; ++i )
    //     std::cout << mP[ i*mNyp2 + 1 ] << ' ';
    // std::cout << std::endl;

    // do
    for ( int iter = 0; iter < mIterMax; ++iter )
    {
        // // Swap the new and old pressures
        // std::swap( mP, mPOld );

        // // Initialize the residual to zero
        // residualSq = 0.;

        // Compute the new value of the pressure 
        for ( int i = 1; i < mNx + 1; ++i )
        {
            for ( int j = 1; j < mNy + 1; ++j )
            {
                unsigned int ij   = i * mNyp2 + j;
                unsigned int ijForRHS = (i-1) * mNy + j-1;
                unsigned int im1j = ( i - 1 ) * mNyp2 + j;
                unsigned int ip1j = ( i + 1 ) * mNyp2 + j;
                unsigned int ijm1 = i * mNyp2 + j - 1;
                unsigned int ijp1 = i * mNyp2 + j + 1;
                // Compute the pressure
                mP[ ij ] = ( 1. - mOmega ) * mP[ ij ] 
                           + 0.25 * mOmega * ( mP[ ip1j ] + mP[ im1j ] +
                                               mP[ ijp1 ] + mP[ ijm1 ] - 
                                               mHSq * mRHS[ ijForRHS ] );

                // mPOld[ ij ] = mP[ ij ];

                // std::cout << mP[ ij ] << ' ' << mP[ ip1j ] << " /// ";
                // std::cout << mPOld[ ij ] << ' ' << mPOld[ ip1j ] << " /// ";

                // Set the boundary values equal to the adjacent ones
                if ( i == 1 )
                    mP[ im1j ] = mP[ ij ];
                else if ( i == mNx )
                    mP[ ip1j ] = mP[ ij ];
                if ( j == 1 )
                    mP[ ijm1 ] = mP[ ij ];
                else if ( j == mNy )
                    mP[ ijp1 ] = mP[ ij ];

                // // Add to the residual
                // double thisRes = mHSqInv * ( mP[ ip1j ] - 2.*mP[ ij ] + mP[ im1j ] + 
                //                           mP[ ijp1 ] - 2.*mP[ ij ] + mP[ ijm1 ] )
                //                  - mRHS[ ijForRHS ];
                // residualSq += thisRes * thisRes;
            }
        }

        // std::cout << mP[ 15 ] << ' ' << mP[ 16 ] << " /// ";
        // std::cout << mPOld[ 15 ] << ' ' << mPOld[ 16 ] << " /// ";

        // Compute the residual
        residualSq = 0.;
        for ( int i = 1; i < mNx + 1; ++i )
        {
            for ( int j = 1; j < mNy + 1; ++j )
            {
                unsigned int ij   = i * mNyp2 + j;
                unsigned int ijForRHS = (i-1) * mNy + j-1;
                unsigned int im1j = ( i - 1 ) * mNyp2 + j;
                unsigned int ip1j = ( i + 1 ) * mNyp2 + j;
                unsigned int ijm1 = i * mNyp2 + j - 1;
                unsigned int ijp1 = i * mNyp2 + j + 1;
                // Add to the residual
                double thisRes = mHSqInv * ( mP[ ip1j ] - 2.*mP[ ij ] + mP[ im1j ] + 
                                             mP[ ijp1 ] - 2.*mP[ ij ] + mP[ ijm1 ] )
                                 - mRHS[ ijForRHS ];
                residualSq += thisRes * thisRes;
            }
        }

        // std::cout << "Iteration " << iter << " of " << mIterMax << 
        //              ", residual squared: " << residualSq << 
        //              " difference " << residualSq - residualSqOld << std::endl;

        // If the change in the residual is small, break
        if ( std::fabs( residualSq - residualSqOld ) < mThreshold )
            break;

        // Save the value of the residual
        residualSqOld = residualSq;
    }
    // while ( residualSq > mThreshold );

    // std::cout << std::endl;
    // std::cout << "Final pressure: ";
    // for ( int i = 1; i < mNx + 1; ++i )
    //     std::cout << mP[ i*mNyp2 + 1 ] << ' ';
    // std::cout << std::endl;
    // std::cout << std::endl << "----------" << std::endl << std::endl;
}

// Update the velocity
void NavierStokesSolver::updateVelocity()
{
    for ( int i = 1; i < mNx + 1; ++i )
    {
        for ( int j = 1; j < mNy + 1; ++j )
        {
            unsigned int ij   = i * mNyp2 + j;
            unsigned int im1j = ( i - 1 ) * mNyp2 + j;
            unsigned int ijm1 = i * mNyp2 + j - 1;

            mU[ ij ] = mUp[ ij ] - mHInv * ( mP[ ij ] - mP[ im1j ] );
            mV[ ij ] = mVp[ ij ] - mHInv * ( mP[ ij ] - mP[ ijm1 ] );
        }
    }
}

// Initialize the file to save the results
void NavierStokesSolver::initOutputFile( const std::string& outfile )
{

}

// Save the results to the file
void NavierStokesSolver::saveStepToFile( const std::string& outFile )
{
    
}

// Close the file to save the results
void NavierStokesSolver::closeOutputFile( const std::string& outfile )
{

}


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
    mThreshold = 0.00001 * mNx * mNy;
    // Maximum number of iterations for the algorithm
    mIterMax = 1000;

    // The simulation is good by default
    mGoodSimulation = true;

    // Scaling for the velocity
    // This will be equal to the maximum value of the velocity in the simulation 
    mVelocityScale = 0.;

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

// Destructor
NavierStokesSolver::~NavierStokesSolver()
{
    // De-allocate the member arrays
    delete[] mU;
    delete[] mV;
    delete[] mRx1;
    delete[] mRx2;
    delete[] mRy1;
    delete[] mRy2;
    delete[] mUp;
    delete[] mVp;
    delete[] mRHS;
    delete[] mP;
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
            if ( mBoundaryConditions->isFluidCell( i-1, j-1 ) )
            {
                mU[ i*mNyp2 + j ] = u0;
                mV[ i*mNyp2 + j ] = v0;
                mP[ i*mNyp2 + j ] = p0;
            }
        }
    }
}

// ------------------------------------------------------
// Main loop

// Run the simulation
void NavierStokesSolver::solve( std::string outFileName )
{
    // Open the output file and write the information in the header
    initOutputFile( outFileName );

    // Initialize the time counter
    auto t0 {std::chrono::high_resolution_clock::now()};

    // Apply the boundary conditions for the first time
    mBoundaryConditions->applyBCs( mU, mV );

    saveStepToFile();

    int count = 0;
    // Iterate over time
    while ( mTotalTime < mTimeSim )
    {
        // If the simulation is not good, stop
        if ( !mGoodSimulation )
            break;

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
        mBoundaryConditions->applyBCs( mU, mV );

        // If the time elapsed is more than a frame, save the current state
        mThisFrameTime += mDeltaT;
        mTotalTime += mDeltaT;
        if ( mThisFrameTime >= mFrameTime )
        {
            mFrameCount += 1;
            mThisFrameTime -= mFrameTime;

            // Save the results to the file
            saveStepToFile();

            std::cout << "Frame " << mFrameCount << " of " << mTotalFrames << " --- " 
                      << mTotalTime << " seconds.\n";
        }
    }

    // Compute the total time of the simulation
	auto t1 {std::chrono::high_resolution_clock::now()};
	std::chrono::duration<double> duration {t1 - t0};
    std::cout << "\nSimulation done in " << duration.count() << " seconds\n";

    // Write the velocity scale and close the output file 
    writeVelocityScale();
    mOutputFile.close();
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

    // Pressure
    mP = new double [ (mNx + 2) * (mNy + 2) ];

    // Initialize all matrices with zeroes
    for ( int i = 0; i < mNxp2; ++i )
    {
        for ( int j = 0; j < mNyp2; ++j )
        {
            mU[ i*mNyp2 + j ] = 0.;
            mV[ i*mNyp2 + j ] = 0.;
            mUp[ i*mNyp2 + j ] = 0.;
            mVp[ i*mNyp2 + j ] = 0.;
            mP[ i*mNyp2 + j ] = 0.;
        }
    }
    for ( int i = 0; i < mNx; ++i )
    {
        for ( int j = 0; j < mNy; ++j )
        {
            mRx1[ i*mNy + j ] = 0.;
            mRy1[ i*mNy + j ] = 0.;
            mRx2[ i*mNy + j ] = 0.;
            mRy2[ i*mNy + j ] = 0.;
            mRHS[ i*mNy + j ] = 0.;
        }
    }
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
    uMax = std::max( std::fabs( uMax ), 0.0001 );
    vMax = std::max( std::fabs( vMax ), 0.0001 );
    double velMax = std::max( uMax, vMax );

    // Take the minimum of the three values 
    mDeltaT = std::min( mReynolds / 4. * mH, 1. / velMax );
    
    // Check that the time step is smaller than the length of the frame
    mDeltaT = std::min( mDeltaT, mFrameTime );

    // Multiply the time step by tau, and by the step h that I left out in the previous line
    mDeltaT *= mTau * mH;
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
            // mUp[ ij ] = mU[ ij ] + mDeltaT * ( 1.5 * mRx[ ijForR ] - 0.5 * mRxOld[ ijForR ] );
            // mVp[ ij ] = mV[ ij ] + mDeltaT * ( 1.5 * mRy[ ijForR ] - 0.5 * mRyOld[ ijForR ] );
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
        }
    }

    // Set the values of the components up[0, j] and vp[i, 0] so that the derivatives
    // next to these boundaries are computed correctly when evaluating the RHS of
    // the Poisson equation 
    for ( int j = 0; j < mNyp2; ++j )
    {
        mUp[ j ] = mU[ j ];
        mUp[ (mNx+1)*mNyp2 + j ] = mU[ (mNx+1)*mNyp2 + j ];
    }
    // for ( int i = 1; i < mNx + 1; ++i )
    for ( int i = 0; i < mNxp2; ++i )
    {
        mVp[ i*mNyp2 ] = mV[ i*mNyp2 ];
        mVp[ i*mNyp2 + mNy+1 ] = mV[ i*mNyp2 + mNy+1 ];
    }

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
            mRHS[ ijForRHS ] = mHInv * ( mUp[ ij ] - mUp[ im1j ] + mVp[ ij ] - mVp[ ijm1 ] ) / mDeltaT;
        }
    }
}

// Solve the Poisson eq. for the pressure 
void NavierStokesSolver::solvePoisson()
{
    // Initialize the L2-norm of the residual
    double residualSq;

    // Set the old pressure to be the average of the squared pressure
    double pSqAcc = 0.;
    for ( int i = 1; i < mNx; ++i )
        for ( int j = 1; j < mNy; ++j )
            pSqAcc += mP[ i*mNyp2 + j ] * mP[ i*mNyp2 + j ];
    pSqAcc = std::sqrt( pSqAcc / ( mNx * mNy ) );
    if ( pSqAcc < 0.001 ) 
        pSqAcc = 1.;
    for ( int i = 0; i < mNxp2; ++i )
        for ( int j = 0; j < mNyp2; ++j )
            mP[ i*mNyp2 + j ] = pSqAcc;

    // Compute the threshold for the SOR algorith based on this value
    double thisTreshold = mThreshold * pSqAcc * mHSqInv / mDeltaT;

    // std::cout << "Average pressure: " << pSqAcc << std::endl;

    int iter = 0;
    for ( ; iter < mIterMax; ++iter )
    {
        // Red/black SOR
        for ( int redBlack = 0; redBlack <=1; ++redBlack )
        {
            // Compute the new value of the pressure 
            for ( int i = 1; i < mNx + 1; ++i )
            {
                for ( int j = 1; j < mNy + 1; ++j )
                {
                    if ( ( i + j ) % 2 != redBlack )
                        continue;
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

                    // Set the boundary values equal to the adjacent ones
                    if ( i == 1 )
                        mP[ im1j ] = 2. * mP[ ij ] - mP[ip1j];
                    else if ( i == mNx )
                        mP[ ip1j ] = 2. * mP[ ij ] - mP[im1j];
                    if ( j == 1 )
                        mP[ ijm1 ] = 2. * mP[ ij ] - mP[ijp1];
                    else if ( j == mNy )
                        mP[ ijp1 ] = 2. * mP[ ij ] - mP[ijm1];
                }
            }
        }

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
        //              // ", squared: " << std::sqrt(residualSq) << 
        //              " difference " << residualSq - residualSqOld << std::endl;

        // If the change in the residual is small, break
        if ( std::fabs( residualSq ) < thisTreshold )
            break;
    }

    if ( iter == mIterMax )
    {
        std::cout << "Reached the maximum number of iterations in SOR.\n";
        // In this case, the simulation is not good
        mGoodSimulation = false;
    }
}

// Update the velocity
void NavierStokesSolver::updateVelocity()
{
    // std::cout << "\n----------------------\n";
    for ( int i = 1; i < mNx + 1; ++i )
    {
        for ( int j = 1; j < mNy + 1; ++j )
        {
            unsigned int ij   = i * mNyp2 + j;
            unsigned int ip1j = ( i + 1 ) * mNyp2 + j;
            unsigned int ijp1 = i * mNyp2 + j + 1;

            mU[ ij ] = mUp[ ij ] - mHInv * mDeltaT * ( mP[ ip1j ] - mP[ ij ] );
            mV[ ij ] = mVp[ ij ] - mHInv * mDeltaT * ( mP[ ijp1 ] - mP[ ij ] );
        }
    }
}

// Initialize the file to save the results
void NavierStokesSolver::initOutputFile( const std::string& outFileName )
{
    // Open the file
    mOutputFile.open( outFileName.c_str(), std::ios::binary | std::ios::out );

    // Write the header
    mOutputFile << "velocity_scale " << "          " << '\n';
    mOutputFile << "dimensions " << mNx << ' ' << mNy << '\n';
    mOutputFile << "cell_size " << mH << '\n';
    mOutputFile << "reynolds_nr " << mReynolds << '\n';
    mOutputFile << "nr_of_frames " << mTimeSim / mFrameTime << '\n';
    mOutputFile << "fps " << 1. / mFrameTime << '\n';

    // Save a map of the grid 
    mBoundaryConditions->drawDomainMap( mOutputFile );
}

// Save the results to the file
void NavierStokesSolver::saveStepToFile()
{
    // If the simulation is not good, do not write anything
    if ( mGoodSimulation )
    {
        // Blank line
        mOutputFile << '\n';

        // Save the velocities u and v, and the pressure, in contiguous form
        for ( int i = 0; i < mNxp2 * mNyp2; ++i )
        {
            mOutputFile << mU[ i ] << ' ';

            if ( std::fabs(mU[ i ]) > mVelocityScale )
                mVelocityScale = std::fabs(mU[i]);
        }
        mOutputFile << '\n';
        for ( int i = 0; i < mNxp2 * mNyp2; ++i )
        {
            mOutputFile << mV[ i ] << ' ';

            if ( std::fabs(mV[ i ]) > mVelocityScale )
                mVelocityScale = std::fabs(mV[i]);
        }
        mOutputFile << '\n';
        for ( int i = 0; i < mNxp2 * mNyp2; ++i )
            mOutputFile << mP[ i ] << ' ';
        mOutputFile << '\n';
    }
}

// Write the velocity scale at the beginning of the file 
void NavierStokesSolver::writeVelocityScale()
{
    // Invert the velocity scale
    mVelocityScale = 2. / mVelocityScale;

    // Move the pointer to the beginning of the file
    mOutputFile.seekp( 15 );
    mOutputFile << std::scientific << std::setprecision(3) << mVelocityScale;
}

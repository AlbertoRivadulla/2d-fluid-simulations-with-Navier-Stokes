#include "videoOutput.h"


// Create a video from the data in a file
void makeVideoFromDataFile( std::string inFileName, std::string outVideoName, int scale )
{
    // Open the file with the data
    std::ifstream dataFile;
    dataFile.open( inFileName );
    if (!dataFile)
        std::cout << "Could not open the file " << inFileName << "\n";

    std::cout << "\n\nWriting a video to the file " << outVideoName << "... " << std::endl;

    // Parameters to read from the file
    int Nx, Ny;
    double h;
    int nrFrames;
    int fps;
    double velScale;
    // String of data to dump
    std::string trashString;

    // Read the parameters
    dataFile >> trashString >> velScale;    // Scaling factor for the velocity
    dataFile >> trashString >> Nx >> Ny;    // Dimensions of the domain
    dataFile >> trashString >> h;           // Size of the cells
    dataFile >> trashString >> trashString; // Reynolds number
    dataFile >> trashString >> nrFrames;    // Number of frames
    dataFile >> trashString >> fps;         // Frames per second

    // Divide the fps by a factor, to make the video slower
    fps = fps / 6;

    // Read the map of the domain 
    std::vector<unsigned int> domainMap ( Nx * Ny, 1 );
    std::string domainMapString;

    dataFile >> domainMapString;

    for ( int i = 0; i < Nx*Ny; ++i )
    {
        if ( domainMapString[ i ] == '0' )
            domainMap[ i ] = 0;
    }

    // Matrices with the data of the two components of the velocity and the pressure 
    double* u = new double[ (Nx + 2) * (Ny + 2) ];
    double* v = new double[ (Nx + 2) * (Ny + 2) ];
    double* p = new double[ (Nx + 2) * (Ny + 2) ];

    // Define the codec and create VideoWriter object.
    // Define the fps to be equal to 10. Also frame size is passed.
    // Notice that the size of the frame here is the opposite as that given when 
    // creating a cv::Mat
    cv::VideoWriter videoWriter ( outVideoName, cv::VideoWriter::fourcc('m','p','4','v'),
                                  fps, cv::Size( Nx * scale, Ny * scale ) );

    // Vector of pixels with double values
    std::vector<cv::Vec3b> pixels ( Nx * Ny );

    // Create an image with this vector as its data
    // Notice the order of the dimensions
    // Integer values, three channels
    cv::Mat frame ( Ny, Nx, CV_8UC3, &pixels[0] );
    // Create a new cv::Mat for resizing the frame
    cv::Mat frameLarge;

    // Create a window to display the image
    cv::namedWindow("Display Image", cv::WINDOW_NORMAL);

    // Write the frames to the video file 
    for ( int thisFrame = 0; thisFrame < nrFrames; ++thisFrame )
    {
        // Read the data for the current frame
        double value;
        for ( int i = 0; i < ( Nx+2 ) * ( Ny+2 ); ++i )
            dataFile >> u[ i ];
        for ( int i = 0; i < ( Nx+2 ) * ( Ny+2 ); ++i )
            dataFile >> v[ i ];
        for ( int i = 0; i < ( Nx+2 ) * ( Ny+2 ); ++i )
            dataFile >> p[ i ];

        // Write the data to the pixels
        for ( int i = 0; i < Nx; ++i )
        {
            for ( int j = 0; j < Ny; ++j )
            {
                // If this cell corresponds to an obstacle, draw it
                if ( domainMap[ i*Ny + j ] == 0 )
                {
                    pixels[ i + (Ny - j - 1)*Nx ] = cv::Vec3b( 0, 255, 0 );
                }
                else
                {
                    // Modulus of the velocity
                    unsigned int ij = ( i + 1 ) * ( Ny + 2 ) + j + 1;
                    double Rval = 255. * velScale * std::sqrt( v[ ij ] * v[ ij ] + u[ ij ] * u[ ij ] );

                    // Absolute value of the vorticity
                    // unsigned int ij = ( i + 1 ) * ( Ny + 2 ) + j + 1;
                    // unsigned int im1j = ( i ) * ( Ny + 2 ) + j + 1;
                    // unsigned int ijm1 = ( i + 1 ) * ( Ny + 2 ) + j;
                    // double Rval = 255. * velScale * ( u[ ij ] - u[ ijm1 ] + v[ ij ] - v[ im1j ] ) / h;
                    // Rval = std::fabs( Rval );

                    pixels[ i + (Ny - j - 1) * Nx ] = cv::Vec3b( Rval, Rval, Rval );
                }
            }
        }

        // If the scale is different than one, resize the frame
        if ( scale != 1 )
        {
            // Resize setting the relative size 
            cv::resize( frame, frameLarge, cv::Size(), scale, scale );
            // Write the frame to the video
            videoWriter.write( frameLarge );
        }
        else 
            // Otherwise write it with the original size
            videoWriter.write( frame );
    }

    // Release the video recorder
    videoWriter.release();

    // Deallocate the arrays
    delete[] u;
    delete[] v;
    delete[] p;

    std::cout << "Done!\n";
}

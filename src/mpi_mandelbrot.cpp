#include <cstdlib>
#include <cstdio>

#include <mpi.h>

#define WORKTAG 1
#define DIETAG 2

// ======================================================
// ======================================================
struct MandelbrotParams
{

    //! image width in pixel
    unsigned int NX;

    //! image height in pixel
    unsigned int NY;

    //! maximum number of iterations in used mandelbrot compute
    unsigned int MAX_ITERS;

    //! maximum number of colors, for scaling the number of iterations
    unsigned int MAX_COLOR;

    //! global domain, lower left corner x coordinate
    double xmin;

    //! global domain, lower left corner y coordinate
    double ymin;

    //! global domain, upper right corner x coordinate
    double xmax;

    //! global domain, upper right corner y coordinate
    double ymax;

    //! distance between two grid points
    double dx;

    //! distance between two grid points
    double dy;

    //! offset to local domain - x coordinate
    int imin;

    //! offset to local domain - y coordinate
    int jmin;

    //! number of pixels in local sub-domain along x axis
    int delta_i;

    //! number of pixels in local sub-domain along y axis
    int delta_j;

    //! default constructor - local subdomain is equal to the entire domain
    MandelbrotParams(int default_size)
            : NX(default_size),
              NY(default_size),
              MAX_ITERS(6000),
              MAX_COLOR(255),
              xmin(-1.7),
              ymin(-1.2),
              xmax(.5),
              ymax(1.2),
              dx(0.0),
              dy(0.0),
              imin(0),
              jmin(0),
              delta_i(NX),
              delta_j(NY)

    {
        dx=(xmax-xmin)/NX;
        dy=(ymax-ymin)/NY;
    }

    MandelbrotParams(int default_size, int imin_, int jmin_, int delta_i_, int delta_j_)
            : MandelbrotParams(default_size)
    {
        imin = imin_;
        jmin = jmin_;
        delta_i = delta_i_;
        delta_j = delta_j_;
    }

}; // struct MandelbrotParams

// ======================================================
// ======================================================
class MandelbrotSet
{

public:
    MandelbrotSet(MandelbrotParams params):
            m_params(params)
    {
        data = new unsigned char[m_params.delta_i*m_params.delta_j];
    }


    ~MandelbrotSet()
    {
        delete[] data;
    }


    void compute()
    {
        auto& imin = m_params.imin;
        auto& jmin = m_params.jmin;
        auto& delta_i = m_params.delta_i;
        auto& delta_j = m_params.delta_j;
        auto& NX = m_params.NX;

        for (int j = jmin; j<jmin+delta_j; ++j)
            for (int i = imin; i<imin+delta_i; ++i)
            {
                // local coordinates
                int il = i - imin;
                int jl = j - jmin;
                data[il + delta_i*jl] = compute_pixel(i, j);
            }
    }

    unsigned char* data;

private:
    MandelbrotParams m_params;


    unsigned char compute_pixel(int Px, int Py) const
    {

        auto& dx = m_params.dx;
        auto& dy = m_params.dy;
        auto& xmin = m_params.xmin;
        auto& ymin = m_params.ymin;

        auto& MAX_ITERS = m_params.MAX_ITERS;
        auto& MAX_COLOR = m_params.MAX_COLOR;

        float x0=xmin+Px*dx;
        float y0=ymin+Py*dy;
        float x=0.0;
        float y=0.0;
        int i = 0;

        for(i=0; x*x+y*y<4.0 and i<MAX_ITERS; i++)
        {
            float xtemp=x*x-y*y+x0;
            y=2*x*y+y0;
            x=xtemp;
        }
        return (unsigned char) MAX_COLOR * i / MAX_ITERS;
    }

}; // class MandelbrotSet

// ======================================================
// ======================================================
void write_screen(unsigned char*            data,
                  const MandelbrotParams&   params)
{
    // print aesthetically, dont read this part
    int xmax=80;
    int ymax=60;
    for(int y=0;y<ymax;y++)
    {
        auto j = y*params.NY/ymax;
        if (j>=params.jmin and j<params.jmin+params.delta_j)
        {
            printf("\n");
            for(int x=0;x<xmax;x++)
            {
                auto i = x*params.NX/xmax;

                // local coordinates
                int il = i - params.imin;
                int jl = j - params.jmin;

                int val = data[il + params.NX*jl];

                if (val==200) printf("&");
                else if (val==42) printf("X");
                else if(val>64) printf("#");
                else if(val>32) printf(":");
                else if(val>8) printf(".");
                else printf(" ");
            }
        }
    }

    printf("\n");

} // write_screen

// ======================================================
// ======================================================

void write_ppm(unsigned char*            data,
               const std::string&        filename,
               const MandelbrotParams&   params)
{
    auto& NX = params.NX;
    auto& NY = params.NY;

    FILE* myfile = fopen(filename.c_str(),"w");

    fprintf(myfile, "P6 %d %d 255\n", NX , NY);
    for(unsigned int i=0; i<NX; ++i)
    {
        for(unsigned int j=0; j<NY; ++j)
        {
            unsigned char pix;
            // create an arbitrary RBG code mapping values taken by imageHost
            if (data[i+NX*j] == 255)
            {
                pix = 0;
                fwrite(&pix,1,1,myfile);
                pix = 0;
                fwrite(&pix,1,1,myfile);
                pix = 0;
                fwrite(&pix,1,1,myfile);
            }
            else if (data[i+NX*j] == 0)
            {
                pix = 0;
                fwrite(&pix,1,1,myfile);
                pix = 0;
                fwrite(&pix,1,1,myfile);
                pix = 55;
                fwrite(&pix,1,1,myfile);
            }
            else
            {
                pix = data[i+NX*j] % 8 * 30;
                fwrite(&pix,1,1,myfile);
                pix = data[i+NX*j] % 8 * 30;
                fwrite(&pix,1,1,myfile);
                pix = 255;
                fwrite(&pix,1,1,myfile);
            }

        }
    }

    fclose(myfile);

} // write_ppm

void master (MandelbrotParams &params)
{
    int nbTask;
    MPI_Status status;

    MPI_Comm_size(MPI_COMM_WORLD, &nbTask);

    // ======================== INIT ======================
    int pos_i = 0;
    int pos_j = 0;

    // First send
    int paramsRank[(nbTask-1)*2];
    int* pos = new int[2];
    int rank;
    int realNbTask = 1;

    for (rank = 1; rank < nbTask; ++rank)
    {
        paramsRank[2*rank-2] = pos_i;
        paramsRank[2*rank-1] = pos_j;

        pos[0] = pos_i;
        pos[1] = pos_j;

        MPI_Send(pos, 2, MPI_INT, rank, WORKTAG, MPI_COMM_WORLD);

        pos_i = pos_i + params.delta_i;

        realNbTask = realNbTask + 1;

        if (pos_i >= params.NX)
        {
            pos_i = 0;
            pos_j = pos_j + params.delta_i;

            if (pos_j >= params.NY)
                break;
        }
    }

    // Kill Non used task
    // Waste of time
    /*for (rank = realNbTask; rank < nbTask; ++rank) {
        MPI_Send(0, 0, MPI_INT, rank, DIETAG, MPI_COMM_WORLD);
    }*/



    //=============== assemble pieces from other MPI processes =======================

    unsigned int NX = params.NX;
    unsigned int NY = params.NY;
    int delta_i = params.delta_i;
    int delta_j = params.delta_j;

    // allocate buffer
    int recv_size = delta_i*delta_j;
    unsigned char * buffer = new unsigned char[recv_size];

    unsigned char *image = new unsigned char[NX*NY];

    while (pos_j < params.NY)
    {
        MPI_Recv(buffer, recv_size, MPI_UNSIGNED_CHAR, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

        rank = status.MPI_SOURCE;

        int imin = paramsRank[2*rank-2];
        int imax = imin+params.delta_i;

        int jmin = paramsRank[2*rank-1];
        int jmax = jmin+params.delta_j;

        // write buffer in global image
        for (int j = jmin; j<jmax; ++j)
        {
            int jl = j - jmin;
            for (int i = imin; i<imax; ++i)
            {
                // local coordinates
                int il = i - imin;
                image[i+NX*j] = buffer[il+params.delta_i*jl];
            }
        }

        paramsRank[2*rank-2] = pos_i;
        paramsRank[2*rank-1] = pos_j;


        pos[0] = pos_i;
        pos[1] = pos_j;

        MPI_Send(pos, 2, MPI_INT, rank, WORKTAG, MPI_COMM_WORLD);

        pos_i = pos_i + params.delta_i;

        if (pos_i >= params.NX)
        {
            pos_i = 0;
            pos_j = pos_j + params.delta_j;
        }
    }
    /*
     * Receive results for pending work requests.
     */
    for (rank = 1; rank < realNbTask; ++rank)
    {
        MPI_Recv(buffer, recv_size, MPI_UNSIGNED_CHAR, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

        int imin = paramsRank[2*status.MPI_SOURCE-2];
        int imax = imin+params.delta_i;

        int jmin = paramsRank[2*status.MPI_SOURCE-1];
        int jmax = jmin+params.delta_j;

        // write buffer in global image
        for (int j = jmin; j<jmax; ++j)
        {
            int jl = j - jmin;
            for (int i = imin; i<imax; ++i)
            {
                // local coordinates
                int il = i - imin;
                image[i+NX*j] = buffer[il+params.delta_i*jl];
            }
        }
    }

    // ======================== FINALIZE ======================
    /*
     * Tell all the workers to exit.
     */
    for (rank = 1; rank < realNbTask; ++rank) {
        MPI_Send(0, 0, MPI_INT, rank, DIETAG, MPI_COMM_WORLD);
    }

    // finally write complete image
    write_ppm(image, "mandelbrot.ppm", params);

    delete[] image;
    delete[] buffer;
    delete[] pos;
}

void worker (MandelbrotParams &params, int myRank)
{
    MPI_Status status;
    int* pos = new int[2];
    
    while (true)
    {
        MPI_Recv(pos, 2, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

        /* Check the tag of the received message */
        if (status.MPI_TAG == DIETAG)
        {
            delete[] pos;
            return;
        }

        params.imin = pos[0];
        params.jmin = pos[1];

        MandelbrotSet mset = MandelbrotSet(params);
        mset.compute();

        // send data to process 0
        int delta_i = params.delta_i;
        int delta_j = params.delta_j;
        int send_size = delta_i*delta_j;

        MPI_Send(mset.data, send_size, MPI_UNSIGNED_CHAR, 0, myRank, MPI_COMM_WORLD);
    }

}

// ======================================================
// ======================================================
int main(int argc, char* argv[])
{

    int global_size = 512;
    int global_cell_size = 32;

    if (argc == 3) {
        global_size = atoi(argv[1]);
        global_cell_size = atoi(argv[2]);
    }

    if (global_cell_size > global_size)
        global_cell_size = global_size;

    if (global_size % global_cell_size != 0)
        return EXIT_FAILURE;

    int cell_size_i = global_cell_size;
    int cell_size_j = global_cell_size;
    int myRank;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

    MandelbrotParams params = MandelbrotParams(global_size,
                              0,
                              0,
                              cell_size_i,
                              cell_size_j);

    if (myRank == 0)
    {
        master(params);
    }
    else
    {
        worker(params, myRank);
    }

    MPI_Finalize();

    return EXIT_SUCCESS;
}

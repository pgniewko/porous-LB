#include "palabos3D.h"
#include "palabos3D.hh"
#include <vector>
#include <cmath>

using namespace plb;
using namespace std;

typedef double T;

#ifdef MRT
    #define DESCRIPTOR descriptors::MRTD3Q19Descriptor
#else
    #define DESCRIPTOR descriptors::D3Q19Descriptor
#endif

// This function object returns a zero velocity, and a pressure which decreases
//   linearly in x-direction. It is used to initialize the particle populations.
class PressureGradient 
{
    public:
        PressureGradient(T deltaP_, plint n_ax_, plint ax_) : deltaP(deltaP_), n_ax(n_ax_), ax(ax_) {}

        void operator() (plint iX, plint iY, plint iZ, T& density, Array<T,3>& velocity) const
        {
            velocity.resetToZero();

            if (ax == 0)
            {
                density = 1. - deltaP * DESCRIPTOR<T>::invCs2 / (T)(n_ax - 1) * (T)iX;
            }
            else if (ax == 1)
            {
                density = 1. - deltaP * DESCRIPTOR<T>::invCs2 / (T)(n_ax - 1) * (T)iY;
            }
            else if (ax == 2)
            {
                density = 1. - deltaP * DESCRIPTOR<T>::invCs2 / (T)(n_ax - 1) * (T)iZ;
            }
        }

    private:
        T deltaP;
        plint n_ax;
        plint ax;
};

Box3D setInlet(plint nx_, plint ny_, plint nz_, plint ax_)
{

    if (ax_ == 0)
    {
        Box3D inlet (0, 0, 0, ny_-1, 0, nz_-1);
        return inlet;
    }
    else if (ax_ == 1)
    {
        Box3D inlet (0, nx_-1, 0, 0, 0, nz_-1);
        return inlet;
    }
    else if (ax_ == 2)
    {
        Box3D inlet (0, nx_-1, 0, ny_-1, 0, 0);
        return inlet;
    }
    else
    {
        Box3D inlet (0, 0, 0, ny_-1, 0, nz_-1);
        return inlet;
    }
}

Box3D setOutlet(plint nx_, plint ny_, plint nz_, plint ax_)
{

    if (ax_ == 0)
    {
        Box3D outlet (nx_-1, nx_-1, 0, ny_-1, 0, nz_-1);
        return outlet;
    }
    else if (ax_ == 1)
    {
        Box3D outlet(0, nx_-1, ny_-1, ny_-1, 0, nz_-1);
        return outlet;
    }
    else if (ax_ == 2)
    {
        Box3D outlet(0, nx_-1, 0, ny_-1, nz_-1, nz_-1);
        return outlet;
    }
    else
    {
        Box3D outlet (nx_-1, nx_-1, 0, ny_-1, 0, nz_-1);
        return outlet;
    }
}

void porousMediaSetup( MultiBlockLattice3D<T,DESCRIPTOR>& lattice,
                       OnLatticeBoundaryCondition3D<T,DESCRIPTOR>* boundaryCondition,
                       MultiScalarField3D<int>& geometry, T deltaP, plint ax)
{
        const plint nx = lattice.getNx();
        const plint ny = lattice.getNy();
        const plint nz = lattice.getNz();

        plint n_ax;
        if (ax == 0)
        {
            n_ax = nx;
        }
        else if (ax == 1)
        {
            n_ax = ny;
        }
        else if (ax == 2)
        {
            n_ax = nz;
        }


        pcerr << "Definition of inlet/outlet." << endl;

        Box3D inlet = setInlet(nx, ny, nz, ax);
        if (ax == 0)
        {
            boundaryCondition->addPressureBoundary0N(inlet, lattice);
        }
        else if (ax == 1)
        {
            boundaryCondition->addPressureBoundary1N(inlet, lattice);
        }
        else if (ax == 2)
        {
            boundaryCondition->addPressureBoundary2N(inlet, lattice);
        }
        setBoundaryDensity(lattice, inlet, (T) 1.);

        Box3D outlet = setOutlet(nx, ny, nz, ax); //(nx-1,nx-1, 1,ny-2, 1,nz-2);
        if (ax == 0)
        {
            boundaryCondition->addPressureBoundary0P(outlet, lattice);
        }
        else if (ax == 1)
        {
            boundaryCondition->addPressureBoundary1P(outlet, lattice);
        }
        else if (ax == 2)
        {
            boundaryCondition->addPressureBoundary2P(outlet, lattice);
        }
        setBoundaryDensity(lattice, outlet, (T) 1. - deltaP*DESCRIPTOR<T>::invCs2);

        pcerr << "Definition of the geometry." << endl;
        // Where "geometry" evaluates to 1, use bounce-back.
        defineDynamics(lattice, geometry, new BounceBack<T,DESCRIPTOR>(), 1);
        // Where "geometry" evaluates to 2, use no-dynamics (which does nothing).
        defineDynamics(lattice, geometry, new NoDynamics<T,DESCRIPTOR>(), 2);

        pcerr << "Initilization of rho and u." << endl;
        initializeAtEquilibrium( lattice, lattice.getBoundingBox(), PressureGradient(deltaP, n_ax, ax) );

        lattice.initialize();
        delete boundaryCondition;
}

void writeVTK(MultiBlockLattice3D<T,DESCRIPTOR>& lattice, plint iter, std::string fout)
{
    VtkImageOutput3D<T> vtkOut(createFileName("vtk_"+fout+"_", iter, 6), 1.);
    vtkOut.writeData<float>(*computeVelocityNorm(lattice), "velocityNorm", 1.);
    vtkOut.writeData<3,float>(*computeVelocity(lattice), "velocity", 1.);

    FileName fn = createFileName("u_3Dfield_"+fout+"_", iter, 6);
    plb_ofstream ofile( ( global::directories().getVtkOutDir() + fn.get() + ".dat").c_str()  );
    Box3D domain = lattice.getBoundingBox();
    std::auto_ptr< MultiScalarField3D< T > > UX = computeVelocityComponent (lattice, domain, 0 );
    std::auto_ptr< MultiScalarField3D< T > > UY = computeVelocityComponent (lattice, domain, 1 );
    std::auto_ptr< MultiScalarField3D< T > > UZ = computeVelocityComponent (lattice, domain, 2 );
    std::auto_ptr< MultiScalarField3D< T > > U = computeVelocityNorm (lattice, domain );

    plint nx = lattice.getNx();
    plint ny = lattice.getNy();
    plint nz = lattice.getNz();

    T u, ux, uy, uz;
    for (plint ix = 0; ix < nx; ix++)
    {
        for (plint iy = 0; iy < ny; iy++)
        { 
            for (plint iz = 0; iz < nz; iz++)
            {
                u  = U->get(ix, iy, iz);
                ux = UX->get(ix, iy, iz);
                uy = UY->get(ix, iy, iz);
                uz = UZ->get(ix, iy, iz);

                if (u > 0)
                    ofile << u << " " << ux << " " << uy << " " << uz << endl;  
                {
                }
            }
        }
    }

    return;
}

T computeTortuosity(MultiBlockLattice3D<T,DESCRIPTOR> &lattice, plint flowComponent )
{

    plint fComponent = flowComponent;
    plint nx = lattice.getNx();
    plint ny = lattice.getNy();
    plint nz = lattice.getNz();
    Box3D pm(0, nx-1, 0, ny-1, 0, nz-1);

    T absvelsum = computeSum(*computeVelocityNorm(lattice, pm));
    T ax_velsum = computeSum(*computeVelocityComponent(lattice, pm, fComponent));

    T t = absvelsum / ax_velsum;
    return t;
}


T computePermeability ( MultiBlockLattice3D<T,DESCRIPTOR>& lattice, T nu, T deltaP, Box3D domain, plint flowComponent )
{
        pcerr << "Computing the permeability." << endl;

        // Compute only the x-direction of the velocity (direction of the flow).
        plint fComponent = flowComponent;
        plint n_ax;
        if (fComponent == 0)
        {
            n_ax = lattice.getNx();
        }
        else if (fComponent == 1)
        {
            n_ax = lattice.getNy();
        }
        else if (fComponent == 2)
        {
            n_ax = lattice.getNz();
        }


        T meanU = computeAverage ( *computeVelocityComponent (lattice, domain, fComponent ) );


        T t = computeTortuosity(lattice, fComponent); 

//        pcout << "Average velocity     = " << meanU                           << endl;
//        pcout << "Lattice viscosity nu = " << nu                              << endl;
//        pcout << "Grad P               = " << deltaP/(T)(n_ax-1)              << endl;
//        pcout << "Permeability         = " << nu*meanU / (deltaP/(T)(n_ax-1)) << endl;
//        pcout << "Tortuosity           = " << t                               << endl;
        
        pcout << meanU                           << " "; // average velocity
        pcout << nu                              << " "; // Lattice viscosity nu
        pcout << deltaP/(T)(n_ax-1)              << " "; // Grad P
        pcout << nu*meanU / (deltaP/(T)(n_ax-1)) << " "; // Permeability
        pcout << t                               << " "; // Tortuosity

        return meanU;
}

int main(int argc, char **argv)
{
        plbInit(&argc, &argv);

        if (argc != 11)
        {
                pcerr << "Error missing some input parameter\n";
                pcerr << "The structure is :\n";
                pcerr << "1. Input file name.\n";
                pcerr << "2. Output directory name.\n";
                pcerr << "4. number of cells in X direction.\n";
                pcerr << "5. number of cells in Y direction.\n";
                pcerr << "6. number of cells in Z direction.\n";
                pcerr << "7. Delta P.\n";
                pcerr << "8. Direction.\n";
                pcerr << "9. Periodicity flag.\n";
                pcerr << "10. Refinemenet level - needed for the max iterations.\n";
                exit (EXIT_FAILURE);
        }
        std::string fNameIn  = argv[1];
        std::string dNameOut = argv[2];
        std::string fNameOut = argv[3];

        const plint nx    = atoi(argv[4]);
        const plint ny    = atoi(argv[5]);
        const plint nz    = atoi(argv[6]);
        const T deltaP    = atof(argv[7]);
        const plint ax    = atoi(argv[8]);
        const plint pbc_  = atoi(argv[9]); 
        const plint ref_k = atoi(argv[10]); 
         
        global::directories().setOutputDir(dNameOut + "/");


        const T omega = 1.0;
        const T nu    = ((T)1/omega-0.5)/DESCRIPTOR<T>::invCs2;
        pcerr << "Cs: " << sqrt( 1.0 / DESCRIPTOR<T>::invCs2 ) << endl;
        pcerr << "Creation of the lattice." << endl;
        

#ifdef MRT
        MultiBlockLattice3D<T,DESCRIPTOR> lattice(nx, ny, nz, new MRTdynamics<T,DESCRIPTOR> ( omega ) );
        pcerr << "MRT dynamics" << endl;
#else
        MultiBlockLattice3D<T,DESCRIPTOR> lattice(nx, ny, nz, new BGKdynamics<T,DESCRIPTOR> ( omega ) );
        pcerr << "BGK dynamics" << endl;
#endif

        // SET PERIODICITY
        if (pbc_ == 1)
        {
            if (ax == 0)
            {
                lattice.periodicity().toggle(0, false);
                lattice.periodicity().toggle(1, true);
                lattice.periodicity().toggle(2, true);
            }
            else if (ax == 1)
            {
                lattice.periodicity().toggle(0, true);
                lattice.periodicity().toggle(1, false);
                lattice.periodicity().toggle(2, true);
            }
            else if (ax == 2)
            {
                lattice.periodicity().toggle(0, true);
                lattice.periodicity().toggle(1, true);
                lattice.periodicity().toggle(2, false);
            }
        }
        else
        {
            // Switch off periodicity.
            lattice.periodicity().toggleAll(false);
        }


        MultiScalarField3D<int> geometry(nx, ny, nz);

        pcerr << "Reading the geometry file." << endl;
        plb_ifstream geometryFile(fNameIn.c_str());
        if (!geometryFile.is_open()) 
        {
                pcerr << "Error: could not open geometry file " << fNameIn << endl;
                return -1;
        }
        geometryFile >> geometry;

        pcerr << "nu = " << nu << endl;
        pcerr << "deltaP = " << deltaP << endl;
        pcerr << "omega = " << omega << endl;
        pcerr << "nx = " << lattice.getNx() << endl;
        pcerr << "ny = " << lattice.getNy() << endl;
        pcerr << "nz = " << lattice.getNz() << endl;

        porousMediaSetup(lattice, createLocalBoundaryCondition3D<T,DESCRIPTOR>(), geometry, deltaP, ax);

        // The value-tracer is used to stop the simulation once is has converged.
        // 1st parameter:velocity
        // 2nd parameter:size
        // 3rd parameters:threshold
        // 1st and second parameters ae used for the length of the time average (size/velocity)
        util::ValueTracer<T> converge(1.0, 1000.0, 1.0e-4);

        pcerr << "Simulation begins" << endl;
        plint iT = 0;

        plint maxT_;

        if (ref_k == 1)
        {
            maxT_ = 100000;
        }
        else if (ref_k == 2)
        {
            maxT_ = 25000;
        }
        
        else if (ref_k == 3)
        {
            maxT_ = 10000;
        }

        const plint maxT = maxT_;

        while(true)
        {
//            if (iT % 200 == 0) 
//            {
//                pcerr << "Iteration " << iT << endl;
//            }
//                
//            if ( (iT+1) % 5000 == 0) 
//            {
//                writeVTK(lattice, iT, fNameOut);
//            }

            lattice.collideAndStream();
            converge.takeValue(getStoredAverageEnergy(lattice), true);

            if (converge.hasConverged())
            {
                break;
            }

            if (iT >= maxT)
            {
                break;
            }
            iT++;
        }

        pcerr << "End of simulation at iteration " << iT << endl;

        pcerr << "Permeability:" << endl << endl;
        computePermeability(lattice, nu, deltaP, lattice.getBoundingBox(), ax);
        pcerr << endl;

        pcerr << "Writing VTK file ..." << endl << endl;
        writeVTK(lattice, iT, fNameOut);
        pcerr << "Finished!" << endl << endl;
}

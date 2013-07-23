/*---Path to ROQ data files---*/
/**---- Rule1: data for Tc=0 ----**/
char p_sizes[]    = "BurstRulesROQ/Rule1/paraSize.txt";
char p_EIM[]      = "BurstRulesROQ/Rule1/DEIM_Indx.txt";
char p_freq_pts[] = "BurstRulesROQ/Rule1/freq_points.txt";
char p_RB_I[]     = "BurstRulesROQ/Rule1/RB_imag.txt";
char p_RB_R[]     = "BurstRulesROQ/Rule1/RB_real.txt";


/*--- Global variables ---*/
gsl_vector *eval;
gsl_matrix *evec;

bool GaussPropSetup= true;
//bool useGaussianProp=false;
bool useGaussianProp=true;

bool SampleAmp=false;
//bool SampleAmp=true;

//bool include_tc = false;
bool include_tc = true;
double fixedtc=0.;

//bool Sample_Tc =false;
bool Sample_Tc =true;

//bool data_save=true;
bool data_save=false;
bool save_lkh_time =false;

int d= 2;       // Initial space-parameter dimension
double snr;     // Signal-to_noise

/*---Information about data stream (for consistency check with ROQ info)---*/
int SampleRate = 20;                        // Should be a power of 2 (in Hz). pow(2,12)
int ObsTime = 45;                           // In seconds
int TDlen = SampleRate*ObsTime;             // Number of time samples
int FDlen = TDlen/2 + 1;                    // Number of non-negative freq samples (+1 for DC)
int SIGlen = 2*FDlen;                       // Number of real+Im signal samples
double fminimum = 0;                        // Smallest non-negative frequency
double fmaximum = (double) SampleRate/2;    // Frequencies resolved up to Nyquist
double deltaF = (double) 1/ObsTime;         // Spacing between frequency samples
double deltaT = (double) 1/SampleRate;
double SigNormSqrd;                         // norm of sig squared, needed for ROQ
char Str_run_number[100];

/* ---Identification with legacy variables--- */
int N = TDlen;
double delta = deltaT;


/*=========== MAIN ===========*/
    // THIS MIGHT BE USEFUL (this was used for creating a directory)
    //system(shell_command);

    // Time or CPU time?
    // http://en.cppreference.com/w/cpp/chrono/c/clock
    //start=clock();

    // Have a switch for ROQ or standard computations

    // initialization of run Parameters

    // Calc C^{-1}

    // Signal amplitude and basis...

    // Use some proposal functions

    //search_start=clock();

    // Do MCMC

    //end=clock();

    // Data save and output statistics

/*---- COMPUTE LIKELIHOODS---*/
// ComputeLogLikelihood_ROQ

// ComputeLogLikelihood

/*==== BEGINNING MCMC ALGORITHM ====*/
void DrawFromPrior(gsl_rng *rng, double *thprop, double *parmin, double *parmax) {
    // A random number
    double ran;

    // WARNING not sure about the amplitude's prior being linear too

    // For each parameter space dimension generate a random number, i.e. select a point
    // randomly in the multi-D space

    // If 
    if ((Sample_Tc) || (SampleAmp)) {
        ran=gsl_rng_uniform(rng);
        thprop[d]=(1.-ran)*parmin[d]+ran*parmax[d];
    } 
    if((Sample_Tc)&&(SampleAmp))
    {
        ran=gsl_rng_uniform(rng);
        thprop[d+1]=(1.-ran)*parmin[d+1]+ran*parmax[d+1];
    }
}

// Calc Fisher matrix

// Use the above matrix to initiate the proposal.
void SetUpGaussianProposal(gsl_matrix *gam)
{
    int Np;
    Np=d;

    if (SampleAmp)
        Np+=1;
    if (Sample_Tc)
        Np+=1;

    eval = gsl_vector_alloc (Np);
    evec = gsl_matrix_alloc (Np, Np);
    gsl_eigen_symmv_workspace * w = gsl_eigen_symmv_alloc (Np);
    gsl_eigen_symmv (gam, eval, evec, w);
    gsl_eigen_symmv_free (w);
    gsl_eigen_symmv_sort (eval, evec, GSL_EIGEN_SORT_ABS_ASC);

    for (int i=0;i<Np;i++)
    {
        // Check what kind of Fisher Matrix was selected
        if (gsl_vector_get(eval,i) < 0.)
        {
            std::cout << "Error - non positive-definite Fisher Matrix detected!" << std::endl;
            exit(0);
        }
        gsl_vector_set(eval,i,1./sqrt(gsl_vector_get(eval,i)));
    }
    GaussPropSetup=true;
}

// Gaussian proposal
// This needs to be initialised
void GaussianProposal(gsl_rng *rng, double *thnow, double *thprop, double *parmin, double *parmax)
{
    double dx;
    if (!GaussPropSetup)
    {
        std::cerr << "Error! Gaussian proposal not initialized." << std::endl;
        exit(0);
    }

    int Np=d;

    if (SampleAmp)
        Np+=1;
    if (Sample_Tc)
        Np+=1;

    bool notfound=true;
    while (notfound)
    {
        for (int i=0;i<Np;i++)
            thprop[i]=thnow[i];

        for (int i=0;i<Np;i++)
        {
            dx=gsl_ran_gaussian(rng,gsl_vector_get(eval,i));

            for (int j=0;j<Np;j++)
                thprop[j]+=dx*gsl_matrix_get(evec,j,i);
        }
        notfound=false;
        for (int i=0;i<Np;i++)
        {
            if ((parmax[i]-thprop[i])*(thprop[i]-parmin[i]) < 0.)
                notfound=true;
        }
    }
    /* for (int i=0;i<Np;i++)
       std::cout << thprop[i] << " ";
       std::cout << std::endl; */
}

/*==== End gaussian proposal ====*/


// MCMC proposal
    // Use either gaussian proposal or draw from prior

// MCMC search
    // Have a switch for ROQ or not ROQ, but maybe this could just be split into several
    // functions?
    //
    // ROQ case
    //  * DrawFromPrior
    //  * Start timing likelihood
    //  * Calculate likelihood
    //  * Stop timing likelihood
    //  * Output some stats
    //
    // Simple case
    //  * DrawFromPrior
    //  * Start timing likelihood
    //  * Calculate likelihood
    //  * Stop timing likelihood
    //  * Output some stats
    //

    // Data saving

    // Metropolis-Hastings
    //
    // for some number of trials, do:
    //  * Propose
    //  * calc likelihood and ratio = exp (likenew - likeold), because likenew is log
    //          likelihood
    //  * accept true
    //  * if ratio < 1. then accept = alpha < ratio, where alpha = rng [0;1]
    //  * If accept, then likeold = likenew, thnow = thprop
    //  * If data save, output to file
/*==== END MCMC ALGORITHM ====*/

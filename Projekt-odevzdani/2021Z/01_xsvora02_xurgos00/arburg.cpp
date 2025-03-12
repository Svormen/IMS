// Implementation from https://c.mql5.com/3/133/Tutorial_on_Burg_smethod_algorithm_recursion.pdf

#include "arburg.h"

// Returns in vector coefficients calculated using Burg algorithm applied to the input source data x
void BurgAlgorithm( std::vector<double> &coeffs, const std::vector<double> &x )
{
    // GET SIZE FROM INPUT VECTORS
    size_t N = x.size() - 1;
    size_t m = coeffs.size();

    // INITIALIZE Ak
    std::vector<double> Ak( m + 1, 0.0 );
    Ak[ 0 ] = 1.0;

    // INITIALIZE f and b
    std::vector<double> f( x );
    std::vector<double> b( x );

    // INITIALIZE Dk
    double Dk = 0.0;
    for ( size_t j = 0; j <= N; j++ )
    {
        Dk += 2.0 * f[ j ] * f[ j ];
    }
    Dk -= f[ 0 ] * f[ 0 ] + b[ N ] * b[ N ];

    // BURG RECURSION
    for ( size_t k = 0; k < m; k++ )
    {
        // COMPUTE MU
        double mu = 0.0;
        for ( size_t n = 0; n <= N - k - 1; n++ )
        {
            mu += f[ n + k + 1 ] * b[ n ];
        }
        mu *= -2.0 / Dk;

        // UPDATE Ak
        for ( size_t n = 0; n <= ( k + 1 ) / 2; n++ )
        {
            double t1 = Ak[ n ] + mu * Ak[ k + 1 - n ];
            double t2 = Ak[ k + 1 - n ] + mu * Ak[ n ];
            Ak[ n ] = t1;
            Ak[ k + 1 - n ] = t2;
        }

        // UPDATE f and b
        for ( size_t n = 0; n <= N - k - 1; n++ )
        {
            double t1 = f[ n + k + 1 ] + mu * b[ n ];
            double t2 = b[ n ] + mu * f[ n + k + 1 ];
            f[ n + k + 1 ] = t1;
            b[ n ] = t2;
        }

        // UPDATE Dk
        Dk = ( 1.0 - mu * mu ) * Dk - f[ k + 1 ] * f[ k + 1 ] - b[ N - k - 1 ] * b[ N - k - 1 ];
    }
    // ASSIGN COEFFICIENTS
    coeffs.assign( ++Ak.begin(), Ak.end() );
}

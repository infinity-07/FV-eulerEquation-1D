#ifndef WENO_RECONSTRUCTOR_HPP
#define WENO_RECONSTRUCTOR_HPP

#include <cmath>
#include <iostream>
#include <iomanip>

inline double weno(double uavemm, double uavem, double uave, double uavep, double uavepp)
{
    // Jiang, G.-S., & Shu, C.-W. (1996).
    // Efficient Implementation of Weighted ENO Schemes.
    // Journal of Computational Physics, 126(1), 202-228. https://doi.org/10.1006/jcph.1996.0130

    // the reconstruction values for each stencil
    const double h0 = (2.0 * uavemm - 7.0 * uavem + 11.0 * uave) / 6.0;
    const double h1 = (-uavem + 5.0 * uave + 2.0 * uavep) / 6.0;
    const double h2 = (2.0 * uave + 5.0 * uavep - uavepp) / 6.0;

    // Linear weights
    const double d0 = 1.0 / 10.0;
    const double d1 = 6.0 / 10.0;
    const double d2 = 3.0 / 10.0;

    // the smoothness indicators
    const double beta0 = 13.0 / 12.0 * pow(uavemm - 2.0 * uavem + uave, 2.0) + 1.0 / 4.0 * pow(uavemm - 4.0 * uavem + 3.0 * uave, 2.0);
    const double beta1 = 13.0 / 12.0 * pow(uavem - 2.0 * uave + uavep, 2.0) + 1.0 / 4.0 * pow(uavem - uavep, 2.0);
    const double beta2 = 13.0 / 12.0 * pow(uave - 2.0 * uavep + uavepp, 2.0) + 1.0 / 4.0 * pow(3.0 * uave - 4.0 * uavep + uavepp, 2.0);

    // Noninear weights
    const double epsilon = 1e-6;
    double alpha0 = d0 / pow((epsilon + beta0), 2.0);
    double alpha1 = d1 / pow((epsilon + beta1), 2.0);
    double alpha2 = d2 / pow((epsilon + beta2), 2.0);

    const double sumAlpha = alpha0 + alpha1 + alpha2;
    const double omega0 = alpha0 / sumAlpha;
    const double omega1 = alpha1 / sumAlpha;
    const double omega2 = alpha2 / sumAlpha;

    // final reconstruction
    const double h = omega0 * h0 + omega1 * h1 + omega2 * h2;

    return h;
}

inline double wenoz(double uavemm, double uavem, double uave, double uavep, double uavepp)
{
    // Borges, R., Carmona, M., Costa, B., & Don, W. S. (2008).
    // An improved weighted essentially non-oscillatory scheme for hyperbolic conservation laws.
    // Journal of Computational Physics, 227(6), 3191-3211. https://doi.org/10.1016/j.jcp.2007.11.038

    // the reconstruction values for each stencil
    const double h0 = (2.0 * uavemm - 7.0 * uavem + 11.0 * uave) / 6.0;
    const double h1 = (-uavem + 5.0 * uave + 2.0 * uavep) / 6.0;
    const double h2 = (2.0 * uave + 5.0 * uavep - uavepp) / 6.0;

    // Linear weights
    const double d0 = 1.0 / 10.0;
    const double d1 = 6.0 / 10.0;
    const double d2 = 3.0 / 10.0;

    // the smoothness indicators
    const double beta0 = 13.0 / 12.0 * pow(uavemm - 2.0 * uavem + uave, 2.0) + 1.0 / 4.0 * pow(uavemm - 4.0 * uavem + 3.0 * uave, 2.0);
    const double beta1 = 13.0 / 12.0 * pow(uavem - 2.0 * uave + uavep, 2.0) + 1.0 / 4.0 * pow(uavem - uavep, 2.0);
    const double beta2 = 13.0 / 12.0 * pow(uave - 2.0 * uavep + uavepp, 2.0) + 1.0 / 4.0 * pow(3.0 * uave - 4.0 * uavep + uavepp, 2.0);

    // Noninear weights
    const double eps = 1e-20;
    const double tau = fabs(beta0 - beta2);
    const double alpha0 = d0 * (1 + tau / (beta0 + eps));
    const double alpha1 = d1 * (1 + tau / (beta1 + eps));
    const double alpha2 = d2 * (1 + tau / (beta2 + eps));

    const double sumAlpha = alpha0 + alpha1 + alpha2;
    const double omega0 = alpha0 / sumAlpha;
    const double omega1 = alpha1 / sumAlpha;
    const double omega2 = alpha2 / sumAlpha;

    // final reconstruction
    const double h = omega0 * h0 + omega1 * h1 + omega2 * h2;

    return h;
}

#endif // WENO_RECONSTRUCTOR_HPP
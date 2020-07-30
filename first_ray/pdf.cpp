#include "pdf.h"

pdf::~pdf() {};

double roughconductor_pdf::pdfVisible(const Vector3f &wi, const Vector3f &m, const hit_record &hrec) const
{
    double cosTheta = dot(hrec.normal, wi);
    if (cosTheta == 0)
        return 0.0;
    return microfacet::smithG1(wi, m, hrec, alphaU, alphaV, distribution_type)
        * std::abs(dot(wi, m)) * microfacet::eval(uvw.fromLocal(m), hrec, alphaU, alphaV, distribution_type) / std::abs(cosTheta);
}
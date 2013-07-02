#include <cmath>
#include <iostream>
#include "vec3.hh"
#include <limits>

#ifndef DBL_EPSILON
#define DBL_EPSILON std::numeric_limits<double>::epsilon()
#endif

namespace cav {  

    // Construction from x y z values:

    Vec3::Vec3 (double x, double y, double z) {
        m_x[0] = x;
        m_x[1] = y;
        m_x[2] = z;
    }

    // Allow empty constructor:

    Vec3::Vec3 () {
        m_x[0] = 0.0;
        m_x[1] = 0.0;
        m_x[2] = 0.0;
    }

    // assignment of one vector to another:

    Vec3& Vec3::operator= (const Vec3& t) {
        m_x[0] = t.x();
        m_x[1] = t.y();
        m_x[2] = t.z();
        return *this;
    }

    // Data access

    // Length of vectors:
    double Vec3::getLengthSquared () const { 
        return m_x[0] * m_x[0] + m_x[1] * m_x[1] + m_x[2] * m_x[2];
    }

    // Output:
    void Vec3::print (std::ostream& os) const {
        os << "(" << m_x[0] << ", " << m_x[1] << ", " << m_x[2] << ")";
    }

    // Unary operators:
    Vec3& Vec3::operator/= (double f) {
        m_x[0] /= f;
        m_x[1] /= f;
        m_x[2] /= f;
        return (*this);
    }

    Vec3& Vec3::operator*= (double f) {
        m_x[0] *= f;
        m_x[1] *= f;
        m_x[2] *= f;
        return (*this);
    }

    Vec3& Vec3::operator+= (const Vec3& v)
    {
        m_x[0] += v.x();
        m_x[1] += v.y();
        m_x[2] += v.z();
        return (*this);
    }

    Vec3& Vec3::operator-= (const Vec3& v)
    {
        m_x[0] -= v.x();
        m_x[1] -= v.y();
        m_x[2] -= v.z();
        return (*this);
    }

    // Boolean operator:
    bool Vec3::operator== (Vec3 x) {
        // all components need to be equal
        // This should be a rigid enough test as it was quite well described at a blog
        // post: http://realtimecollisiondetection.net/blog/?p=89
        bool r = true;
        for (unsigned i = 0; i < 3 and r; i++) {
            double absTol = DBL_EPSILON;
            double relTol = DBL_EPSILON * fmax(fabs(m_x[i]), fabs(x[i]));
            r = r and fabs(m_x[i] - x[i]) <= fmax(absTol, relTol * fmax(fabs(m_x[i]), fabs(x[i])));
        }

        return r;
    }


} // end namespace cav

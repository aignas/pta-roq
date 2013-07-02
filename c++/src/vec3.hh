#ifndef cavlib_vector_hh
#define cavlib_vector_hh

// Define the Vec3 class

// A Vec3 is a vector in three-space. The implementation is all in
// this header file and therefore will be *inlined* for efficiency by
// the compiler (if it can...)
//
// This was edited by me as some on the internet strongly argue, that
// only one line functions should be specified in the class definition
// itself so as to separate the interface and the implementation, which
// always make it clearer what the class does.

#include <cmath>
#include <iostream>

namespace cav {  

class Vec3
{
    private:
        // The only data needed are three double precision values x,y,z:
        double m_x[3];

    public:
        // Construction from x y z values:
        Vec3 (double, double, double);

        // Allow empty constructor:
        Vec3 () ;

        // Destruction: this is not necessary for this class, but we provide
        // an empty destructor anyway as an example
        ~Vec3 () {}

        // assignment of one vector to another:
        Vec3& operator= (const Vec3& t);

        // copy constructor: 
        Vec3 (const Vec3& t) { *this = t; }

        // Data access
  
        // x() y() z() are const, so cannot be used to modify the values:
        double x() const { return m_x[0]; }
        double y() const { return m_x[1]; }
        double z() const { return m_x[2]; }
  
        // ... but one can change values via [] which returns a reference.
        // Note - no bounds checking is done.
        double& operator[] (int i) { return m_x[i]; }
  
        // Length of vectors:
        double getLengthSquared () const;
        double getLength () const {return sqrt(getLengthSquared());}
  
        // Output:
        void print (std::ostream& os) const;
  
        // Normalisation, in place:
        Vec3& normalise () {return *this /= getLength();}
  
        // Unary operators:
        Vec3& operator/= (double f);
        Vec3& operator*= (double f);
        Vec3& operator+= (const Vec3& v);
        Vec3& operator-= (const Vec3& v);
        
        bool operator== (Vec3 x);
};  // end of Vec3 class 

/*
    // I/O of Vec3 objects:
    std::ostream& operator<< (std::ostream& os, const Vec3& a);

    // unary minus negates the vector: 
    Vec3 operator- (const Vec3& v);
    // Addition of vectors:
    Vec3 operator+ (const Vec3& v1, const Vec3& v2);
    // Subtraction:
    Vec3 operator- (const Vec3& v1, const Vec3& v2);
    // Scaling:
    Vec3 operator* (const Vec3& v, double f);
    Vec3 operator* (double f, const Vec3& v);
    Vec3 operator/ (const Vec3& v, double f);
    // inner ("dot") product: v1 % v2
    double operator% (const Vec3& v1, const Vec3& v2);
    // cross product:  v1 * v2
    Vec3 operator* (const Vec3& v1, const Vec3& v2);
    // Return a normalised vector:
    Vec3 normalise (const Vec3& v);
    // Rotate a vector r *clockwise* about given axis by angle phi:
    Vec3 rotate (const Vec3& r, const Vec3& rotation_axis, double phi);
 */


    // I/O of Vec3 objects:
    inline std::ostream& operator<< (std::ostream& os, const Vec3& a) {
        a.print (os);
        return os;
    }

    // Inline non-member function utilities for Vec3 objects:

    // unary minus negates the vector: 
    inline Vec3 operator- (const Vec3& v) {
        return Vec3(-v.x(), -v.y(), -v.z());
    }

    // Addition of vectors:
    inline Vec3 operator+ (const Vec3& v1, const Vec3& v2) {
        return Vec3 (v1.x() + v2.x(), v1.y() + v2.y(), v1.z() + v2.z());
    }

    // Subtraction:
    inline Vec3 operator- (const Vec3& v1, const Vec3& v2 ) {
        return Vec3 (v1.x() - v2.x(), v1.y() - v2.y(), v1.z() - v2.z());
    }

    // Scaling:
    inline Vec3 operator* ( const Vec3& v, double f ) {
        return Vec3( f * v.x(), f * v.y(), f * v.z() );
    }

    inline Vec3 operator* (double f, const Vec3& v ) {
        return Vec3( f * v.x(), f * v.y(), f * v.z() );
    }

    inline Vec3 operator/ (const Vec3& v, double f ) {
        return Vec3( v.x() / f, v.y() / f, v.z() / f );
    }

    // inner ("dot") product: v1 % v2
    inline double operator% (const Vec3& v1, const Vec3& v2 ) {
        return v1.x() * v2.x() + v1.y() * v2.y() + v1.z() * v2.z();
    }

    // cross product:  v1 * v2
    inline Vec3 operator* (const Vec3& v1, const Vec3& v2 ) {
        return Vec3(v1.y() * v2.z() - v1.z() * v2.y(),
                    v1.z() * v2.x() - v1.x() * v2.z(),
                    v1.x() * v2.y() - v1.y() * v2.x());
    }

    // Return a normalised vector:
    inline Vec3 normalise (const Vec3& v) {
        return v / v.getLength();
    }

    // Rotate a vector r *clockwise* about given axis by angle phi:

    inline Vec3 rotate (const Vec3& r, const Vec3& rotation_axis, double phi ) {
        Vec3 n_hat = normalise( rotation_axis );

        return r * cos( phi ) + 
            ( ( (n_hat % r ) * ( 1 - cos(phi) ) ) * n_hat ) -
            ( ( r * n_hat ) * sin(phi) );
    }


} // end namespace cav

#endif


// Copyright (c) 2009-2018 The Regents of the University of Michigan
// This file is part of the HOOMD-blue project, released under the BSD 3-Clause License.

#include "hoomd/HOOMDMath.h"
#include "hoomd/BoxDim.h"
#include "HPMCPrecisionSetup.h"
#include "hoomd/VectorMath.h"
#include "ShapeSphere.h"
#include "XenoCollide3D.h"
#include "ShapeConvexPolyhedron.h"
#include "hoomd/AABB.h"

#ifndef __SHAPE_FACETED_SPHERE_H__
#define __SHAPE_FACETED_SPHERE_H__

/*! \file ShapeFacetedEllipsoid.h
    \brief Defines the faceted ellipsoid shape
*/

/*! The faceted ellispoid is defined by the intersection of an ellipsoid of half axes a,b,c with planes given by the face
 * normals n and offsets b, obeying the equation:
 *
 * r.n + b <= 0
 *
 * Intersections of planes within the ellipsoid volume are accounted for.
 * Internally, the overlap check works by computing the support function for an ellipsoid deformed into a unit sphere.
 *
 */
// need to declare these class methods with __device__ qualifiers when building in nvcc
// DEVICE is __device__ when included in nvcc and blank when included into the host compiler
#ifdef NVCC
#define DEVICE __device__
#define HOSTDEVICE __host__ __device__
#else
#define DEVICE
#define HOSTDEVICE
#endif

namespace hpmc
{

namespace detail
{

//! Data structure for intersection planes
/*! \ingroup hpmc_data_structs */
struct faceted_ellipsoid_params : param_base
    {
    //! Empty constructor
    faceted_ellipsoid_params()
        : a(1.0), b(1.0), c(1.0), N(0), ignore(1)
        { }

    #ifndef NVCC
    faceted_ellipsoid_params(unsigned int n_facet, bool managed )
        : a(1.0), b(1.0), c(1.0), N(n_facet), ignore(0)
        {
        n = ManagedArray<vec3<OverlapReal> >(n_facet, managed);
        offset = ManagedArray<OverlapReal> (n_facet, managed);
        }
    #endif

    poly3d_verts verts;           //!< Vertices of the polyhedron
    poly3d_verts additional_verts;//!< Vertices of the polyhedron edge-sphere intersection
    ManagedArray<vec3<OverlapReal> > n;              //!< Normal vectors of planes
    ManagedArray<OverlapReal> offset;                //!< Offset of every plane
    OverlapReal a;                                   //!< First half-axis
    OverlapReal b;                                   //!< Second half-axis
    OverlapReal c;                                   //!< Third half-axis
    vec3<OverlapReal> origin;                        //!< Origin shift
    unsigned int N;                                  //!< Number of cut planes
    unsigned int ignore;                             //!< Bitwise ignore flag for stats, overlaps. 1 will ignore, 0 will not ignore
                                                     //   First bit is ignore overlaps, Second bit is ignore statistics

    //! Load dynamic data members into shared memory and increase pointer
    /*! \param ptr Pointer to load data to (will be incremented)
        \param available_bytes Size of remaining shared memory allocation
     */
    HOSTDEVICE void load_shared(char *& ptr, unsigned int &available_bytes) const
        {
        n.load_shared(ptr,available_bytes);
        offset.load_shared(ptr,available_bytes);
        verts.load_shared(ptr,available_bytes);
        additional_verts.load_shared(ptr, available_bytes);
        }

    #ifdef ENABLE_CUDA
    //! Attach managed memory to CUDA stream
    void attach_to_stream(cudaStream_t stream) const
        {
        n.attach_to_stream(stream);
        offset.attach_to_stream(stream);
        verts.attach_to_stream(stream);
        additional_verts.attach_to_stream(stream);
        }
    #endif
    } __attribute__((aligned(32)));

//! Support function for ShapeFacetedEllipsoid
/* \ingroup minkowski
*/
class SupportFuncFacetedEllipsoid
    {
    public:
        //! Construct a support function for a faceted sphere
        /*! \param _params Parameters of the faceted sphere
        */
        DEVICE SupportFuncFacetedEllipsoid(const faceted_ellipsoid_params& _params,
            const OverlapReal& _sweep_radius)
            : params(_params), sweep_radius(_sweep_radius)
            {
            }

        DEVICE inline bool isInside(const vec3<OverlapReal>& v,unsigned int plane) const
            {
            // is this vertex masked by a plane?
            vec3<OverlapReal> np = params.n[plane];

            // transform plane normal into the coordinate system of the unit sphere
            np.x *= params.a; np.y *= params.b; np.z *= params.c;
            OverlapReal b = params.offset[plane];

            // is current supporting vertex inside the half-space defined by this plane?
            return (dot(np,v) + b <= OverlapReal(0.0));
            }


        //! Compute the support function
        /*! \param n Normal vector input (in the local frame)
            \returns Local coords of the point furthest in the direction of n
        */
        DEVICE vec3<OverlapReal> operator() (const vec3<OverlapReal>& n_in) const
            {
            // transform support direction into coordinate system of the unit sphere
            vec3<OverlapReal> n(n_in);
            n.x *= params.a; n.y *= params.b; n.z *= params.c;

            OverlapReal nsq = dot(n,n);
            vec3<OverlapReal> n_sphere = n*fast::rsqrt(nsq);

            vec3<OverlapReal> max_vec = n_sphere;
            bool have_vertex = true;

            unsigned int n_planes = params.N;
            for (unsigned int i = 0; i < n_planes; ++i)
                {
                if (! isInside(max_vec,i))
                    {
                    have_vertex = false;
                    break;
                    }
                }

            // iterate over intersecting planes
            for (unsigned int i = 0; i < n_planes; i++)
                {
                vec3<OverlapReal> n_p = params.n[i];
                // transform plane normal into the coordinate system of the unit sphere
                n_p.x *= params.a; n_p.y *= params.b; n_p.z *= params.c;

                OverlapReal np_sq = dot(n_p,n_p);
                OverlapReal b = params.offset[i];

                // compute supporting vertex on intersection boundary (circle)
                // between plane and sphere
                OverlapReal alpha = dot(n_sphere,n_p);
                OverlapReal arg = (OverlapReal(1.0)-alpha*alpha/np_sq);
                vec3<OverlapReal> v;
                if (arg >= OverlapReal(SMALL))
                    {
                    OverlapReal arg2 = OverlapReal(1.0)-b*b/np_sq;
                    OverlapReal invgamma = fast::sqrt(arg2/arg);

                    // Intersection vertex that maximizes support function
                    v = invgamma*(n_sphere-alpha/np_sq*n_p)-n_p*b/np_sq;
                    }
                else
                    {
                    // degenerate case
                    v = -b*n_p/np_sq;
                    }

                bool valid = true;
                for (unsigned int j = 0; j < n_planes; ++j)
                    {
                    if (i!=j && !isInside(v,j))
                        {
                        valid = false;
                        break;
                        }
                    }

                if (valid && (!have_vertex || dot(v,n) > dot(max_vec,n)))
                    {
                    max_vec = v;
                    have_vertex = true;
                    }
                }

            // plane-plane-sphere intersection vertices
            if (params.additional_verts.N)
                {
                detail::SupportFuncConvexPolyhedron s(params.additional_verts);
                vec3<OverlapReal> v = s(n);

                if (! have_vertex || dot(v,n)>dot(max_vec,n))
                    {
                    max_vec = v;
                    have_vertex = true;
                    }
                }

            // plane-plane intersections from user input
            if (params.verts.N)
                {
                detail::SupportFuncConvexPolyhedron s(params.verts);
                vec3<OverlapReal> v = s(n);
                if (dot(v,v) <= OverlapReal(1.0) && (!have_vertex || dot(v,n) > dot(max_vec,n)))
                    {
                    max_vec = v;
                    have_vertex = true;
                    }
                }

            // transform vertex on unit sphere back onto ellipsoid surface
            max_vec.x *= params.a; max_vec.y *= params.b; max_vec.z *= params.c;

            // origin shift
            max_vec -= params.origin;

            // extend out by sweep radius
            return max_vec + (sweep_radius * fast::rsqrt(dot(n_in,n_in))) * n_in;
            }

    private:
        const faceted_ellipsoid_params& params;      //!< Definition of faceted sphere
        const OverlapReal sweep_radius;             //!< The radius of a sphere sweeping the shape
    };



} // end namespace detail


//! Faceted sphere shape template
/*! ShapeFacetedEllipsoid implements IntegratorHPMC's shape protocol for a sphere that is truncated
    by a set of planes, defined through their plane equations n_i*x = n_i^2.

    The parameter defining the sphere is just a single Scalar, the sphere radius.

    \ingroup shape
*/
struct ShapeFacetedEllipsoid
    {
    //! Define the parameter type
    typedef detail::faceted_ellipsoid_params param_type;

    //! Initialize a shape at a given position
    DEVICE ShapeFacetedEllipsoid(const quat<Scalar>& _orientation, const param_type& _params)
        : orientation(_orientation), params(_params)
        { }

    //! Does this shape have an orientation
    DEVICE bool hasOrientation() { return params.N > 0; }

    //!Ignore flag for acceptance statistics
    DEVICE bool ignoreStatistics() const { return params.ignore; }

    //! Get the circumsphere diameter
    DEVICE OverlapReal getCircumsphereDiameter() const
        {
        return OverlapReal(2)*detail::max(params.a, detail::max(params.b, params.c));
        }

    //! Get the in-sphere radius
    DEVICE OverlapReal getInsphereRadius() const
        {
        return 0.0;
        }

    //! Support function of the shape (in local coordinates), used in getAABB
    /*! \param n Vector to query support function (must be normalized)
    */
    DEVICE vec3<Scalar> sfunc(vec3<Scalar> n) const
        {
        vec3<Scalar> numerator(params.a*params.a*n.x, params.b*params.b*n.y, params.c*params.c*n.z);
        vec3<Scalar> dvec(params.a*n.x, params.b*n.y, params.c*n.z);
        return numerator / fast::sqrt(dot(dvec, dvec));
        }

    //! Return the bounding box of the shape in world coordinates
    DEVICE detail::AABB getAABB(const vec3<Scalar>& pos) const
        {
        //OverlapReal max_axis = detail::max(axes.x, detail::max(axes.y, axes.z));
        //return detail::AABB(pos, max_axis);

        // use support function of the ellipsoid to determine the furthest extent in each direction
        vec3<Scalar> e_x(1,0,0);
        vec3<Scalar> e_y(0,1,0);
        vec3<Scalar> e_z(0,0,1);
        vec3<Scalar> s_x = rotate(orientation, sfunc(rotate(conj(orientation),e_x)));
        vec3<Scalar> s_y = rotate(orientation, sfunc(rotate(conj(orientation),e_y)));
        vec3<Scalar> s_z = rotate(orientation, sfunc(rotate(conj(orientation),e_z)));

        // translate out from the position by the furthest extent
        vec3<Scalar> upper(pos.x + s_x.x, pos.y + s_y.y, pos.z + s_z.z);
        // the furthest extent is symmetrical
        vec3<Scalar> lower(pos.x - s_x.x, pos.y - s_y.y, pos.z - s_z.z);

        return detail::AABB(lower, upper);
        }

    //! Return a tight fitting OBB
    DEVICE detail::OBB getOBB(const vec3<Scalar>& pos) const
        {
        // just use the AABB for now
        return detail::OBB(getAABB(pos));
        }

    //! Returns true if this shape splits the overlap check over several threads of a warp using threadIdx.x
    HOSTDEVICE static bool isParallel() { return false; }

    //! Returns true if the overlap check supports sweeping both shapes by a sphere of given radius
    HOSTDEVICE static bool supportsSweepRadius()
        {
        return true;
        }

    /*!
     * Generate the intersections points of polyhedron edges with the sphere
     */
    DEVICE static void initializeVertices(param_type& _params, bool managed)
        {
        #ifndef NVCC
        _params.additional_verts = detail::poly3d_verts(2*_params.N*_params.N, managed);
        _params.additional_verts.diameter = OverlapReal(2.0); // for unit sphere
        _params.additional_verts.N = 0;

        // iterate over unique pairs of planes
        for (unsigned int i = 0; i < _params.N; ++i)
            {
            vec3<OverlapReal> n_p(_params.n[i]);
            // transform plane normal into the coordinate system of the unit sphere
            n_p.x *= _params.a; n_p.y *= _params.b; n_p.z *= _params.c;

            OverlapReal b(_params.offset[i]);

            for (unsigned int j = i+1; j < _params.N; ++j)
                {
                vec3<OverlapReal> np2(_params.n[j]);
                // transform plane normal into the coordinate system of the unit sphere
                np2.x *= _params.a; np2.y *= _params.b; np2.z *= _params.c;

                OverlapReal b2 = _params.offset[j];
                OverlapReal np2_sq = dot(np2,np2);

                // determine intersection line between plane i and plane j

                // point on intersection line closest to center of sphere
                OverlapReal dotp = dot(np2,n_p);
                OverlapReal denom = dotp*dotp-dot(n_p,n_p)*np2_sq;
                OverlapReal lambda0 = OverlapReal(2.0)*(b2*dotp - b*np2_sq)/denom;
                OverlapReal lambda1 = OverlapReal(2.0)*(-b2*dot(n_p,n_p)+b*dotp)/denom;

                vec3<OverlapReal> r = -(lambda0*n_p+lambda1*np2)/OverlapReal(2.0);

                // if the line does not intersect the sphere, ignore
                if (dot(r,r) > OverlapReal(1.0))
                    continue;

                // the line intersects with the sphere at two points, one of
                // maximizes the support function
                vec3<OverlapReal> c01 = cross(n_p,np2);
                OverlapReal s = fast::sqrt((OverlapReal(1.0)-dot(r,r))*dot(c01,c01));

                OverlapReal t0 = (-dot(r,c01)-s)/dot(c01,c01);
                OverlapReal t1 = (-dot(r,c01)+s)/dot(c01,c01);

                vec3<OverlapReal> v1 = r+t0*c01;
                vec3<OverlapReal> v2 = r+t1*c01;

                // check first point
                bool allowed = true;
                for (unsigned int k = 0; k < _params.N; ++k)
                    {
                    if (k == i || k == j) continue;

                    vec3<OverlapReal> np3(_params.n[k]);
                    // transform plane normal into the coordinate system of the unit sphere
                    np3.x *= _params.a; np3.y *= _params.b; np3.z *= _params.c;

                    OverlapReal b3(_params.offset[k]);

                    // is this vertex inside the volume bounded by all halfspaces?
                    if (dot(np3,v1) + b3 > OverlapReal(0.0))
                        {
                        allowed = false;
                        break;
                        }
                    }

                if (allowed && dot(v1,v1) <= OverlapReal(1.0)*(1+SMALL))
                    {
                    _params.additional_verts.x[_params.additional_verts.N] = v1.x;
                    _params.additional_verts.y[_params.additional_verts.N] = v1.y;
                    _params.additional_verts.z[_params.additional_verts.N] = v1.z;
                    _params.additional_verts.N++;
                    }

                // check second point
                allowed = true;
                for (unsigned int k = 0; k < _params.N; ++k)
                    {
                    if (k == i || k == j) continue;

                    vec3<OverlapReal> np3(_params.n[k]);
                    // transform plane normal into the coordinate system of the unit sphere
                    np3.x *= _params.a; np3.y *= _params.b; np3.z *= _params.c;

                    OverlapReal b3(_params.offset[k]);

                    // is this vertex inside the volume bounded by all halfspaces?
                    if (dot(np3,v2) + b3 > OverlapReal(0.0))
                        {
                        allowed = false;
                        break;
                        }
                    }

                if (allowed && dot(v2,v2) <= (OverlapReal(1.0)+SMALL))
                    {
                    _params.additional_verts.x[_params.additional_verts.N] = v2.x;
                    _params.additional_verts.y[_params.additional_verts.N] = v2.y;
                    _params.additional_verts.z[_params.additional_verts.N] = v2.z;
                    _params.additional_verts.N++;
                    }
                }
            }
        #endif
        }

    quat<Scalar> orientation;    //!< Orientation of the sphere (unused)

    const param_type& params;           //!< Faceted sphere parameters
    };

//! Check if circumspheres overlap
/*! \param r_ab Vector defining the position of shape b relative to shape a (r_b - r_a)
    \param a first shape
    \param b second shape
    \returns true if the circumspheres of both shapes overlap

    \ingroup shape
*/
DEVICE inline bool check_circumsphere_overlap(const vec3<Scalar>& r_ab, const ShapeFacetedEllipsoid& a,
    const ShapeFacetedEllipsoid &b)
    {
    OverlapReal DaDb = a.getCircumsphereDiameter() + b.getCircumsphereDiameter();
    vec3<OverlapReal> dr(r_ab);

    return (dot(dr,dr) <= DaDb*DaDb/OverlapReal(4.0));
    }

//! Overlap of faceted spheres
/*! \param r_ab Vector defining the position of shape b relative to shape a (r_b - r_a)
    \param a first shape
    \param b second shape
    \param err in/out variable incremented when error conditions occur in the overlap test
    \param sweep_radius Additional sphere radius to sweep the shapes with
    \returns true when *a* and *b* overlap, and false when they are disjoint

    \ingroup shape
*/
template <>
DEVICE inline bool test_overlap<ShapeFacetedEllipsoid, ShapeFacetedEllipsoid>(const vec3<Scalar>& r_ab, const ShapeFacetedEllipsoid& a, const ShapeFacetedEllipsoid& b,
    unsigned int& err, Scalar sweep_radius_a, Scalar sweep_radius_b)
    {
    vec3<OverlapReal> dr(r_ab);

    // prepend with ellipsoid overlap check?

    OverlapReal DaDb = a.getCircumsphereDiameter() + b.getCircumsphereDiameter();
    return detail::xenocollide_3d(detail::SupportFuncFacetedEllipsoid(a.params, sweep_radius_a),
                           detail::SupportFuncFacetedEllipsoid(b.params, sweep_radius_b),
                           rotate(conj(quat<OverlapReal>(a.orientation)), dr + rotate(quat<OverlapReal>(b.orientation),b.params.origin))-a.params.origin,
                           conj(quat<OverlapReal>(a.orientation))* quat<OverlapReal>(b.orientation),
                           DaDb/2.0,
                           err);

    /*
    return detail::gjke_3d(detail::SupportFuncFacetedEllipsoid(a.params),
                           detail::SupportFuncFacetedEllipsoid(b.params),
                           rotate(conj(quat<OverlapReal>(a.orientation)), dr),
                           conj(quat<OverlapReal>(a.orientation))* quat<OverlapReal>(b.orientation),
                           DaDb/2.0,
                           err);
    */
    }

}; // end namespace hpmc

#endif //__SHAPE_FACETED_SPHERE_H__

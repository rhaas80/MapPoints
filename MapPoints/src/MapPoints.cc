#include <iostream>

#include "cctk.h"
#include "cctk_Parameters.h"

#include "util_ErrorCodes.h"
#include "util_Table.h"

#include "vect.hh"

#include "carpet.hh"


namespace MapPointsNames {

  using namespace Carpet;

  extern "C" CCTK_INT
  MapPoints_MapPoints (CCTK_POINTER_TO_CONST const cctkGH_,
                       CCTK_INT const N_dims,
                       CCTK_INT const param_table_handle,
                       CCTK_INT const coord_system_handle,
                       CCTK_INT const N_interp_points,
                       CCTK_INT const interp_coords_type_code,
                       CCTK_POINTER_TO_CONST const coords_list [],
                       CCTK_POINTER procs,
                       CCTK_POINTER rlev);
  
  static int
  extract_parameter_table_options (cGH const * const cctkGH,
                                   int const param_table_handle,
                                   int const N_interp_points,
                                   bool& have_source_map,
                                   vector<CCTK_INT> & source_map);


  static void
  map_points (cGH const * const cctkGH,
            int const coord_system_handle,
            int const coord_group,
            int const ml,
            int const minrl,
            int const maxrl,
            int const maxncomps,
            int const N_dims,
            int const npoints,
            vector<CCTK_INT> & source_map,
            void const * const coords_list[],
            CCTK_INT procs[],
            CCTK_INT rlev[]);




  extern "C" CCTK_INT
  MapPoints_MapPoints (CCTK_POINTER_TO_CONST const cctkGH_,
                       CCTK_INT const N_dims,
                       CCTK_INT const param_table_handle,
                       CCTK_INT const coord_system_handle,
                       CCTK_INT const N_interp_points,
                       CCTK_INT const interp_coords_type_code,
                       CCTK_POINTER_TO_CONST const coords_list [],
                       CCTK_POINTER procs,
                       CCTK_POINTER rlev)
  {
    DECLARE_CCTK_PARAMETERS;

    static Timer * timer_CDI = NULL;
    if (not timer_CDI) {
      timer_CDI = new Timer ("MapPoints::MapPoints_MapPoints");
    }
    timer_CDI->start();

    cGH const * const cctkGH = static_cast<cGH const *> (cctkGH_);
    assert (cctkGH);
    assert (0 <= N_dims and N_dims <= dim);

    // Check input arrays
    int coord_group = CCTK_GroupIndex ("grid::coordinates"); // use pseudo group
    assert (coord_group >= 0);

    assert (N_interp_points >= 0);
    assert (coords_list);
    for (int d = 0; d < N_dims; ++d) {
      assert (N_interp_points == 0 or coords_list[d]);
    }

    if (interp_coords_type_code != CCTK_VARIABLE_REAL) {
      CCTK_WARN (CCTK_WARN_ABORT, "CarpetInterp does not support interpolation "
                 "coordinates other than datatype CCTK_VARIABLE_REAL");
    }

    // Check output arrays
    assert (N_interp_points == 0 or procs != 0);
    assert (N_interp_points == 0 or rlev != 0);

    // Multiple convergence levels are not supported
    assert (mglevels == 1);
    int const ml = 0;

    //////////////////////////////////////////////////////////////////////
    // Extract parameter table options:
    //   - source map
    //////////////////////////////////////////////////////////////////////
    vector<CCTK_INT> source_map (N_interp_points);
    bool have_source_map;

    {
      int const iret = extract_parameter_table_options
        (cctkGH,
         param_table_handle,
         N_interp_points, 
         have_source_map, 
         source_map);
      if (iret < 0) {
        timer_CDI->stop (0);
        return iret;
      }
    }

    // Find range of refinement levels
    assert (maps > 0);
    for (int m=1; m<maps; ++m) {
      assert (arrdata.AT(coord_group).AT(0).hh->reflevels() ==
              arrdata.AT(coord_group).AT(m).hh->reflevels());
    }
    int const minrl = 0;
    int const maxrl = arrdata.AT(coord_group).AT(0).hh->reflevels();

    // Find maximum number of components over all levels and maps
    int maxncomps = 0;
    for (int rl=minrl; rl<maxrl; ++rl) {
      for (int m=0; m<maps; ++m) {
        maxncomps = max(maxncomps, arrdata.AT(coord_group).AT(m).hh->components(rl));
      }
    }

    // Each point from coord_list is mapped onto the processor
    // that owns it (dstprocs)
    // Also accumulate the number of points per processor (sendcnt)
    map_points (cctkGH, coord_system_handle, coord_group,
                ml, minrl, maxrl, maxncomps, N_dims, N_interp_points,
                source_map,
                coords_list, 
                static_cast<CCTK_INT * const>(procs), 
                static_cast<CCTK_INT * const>(rlev));

    // Done.
    timer_CDI->stop (0);
    return 0;
  }




  static int
  extract_parameter_table_options (cGH const* const cctkGH,
                                   int const param_table_handle,
                                   int const N_interp_points,
                                   bool& have_source_map,
                                   vector<CCTK_INT>& source_map)
  {
    DECLARE_CCTK_PARAMETERS;

    int iret;

    // Find source map
    assert ((int)source_map.size() == N_interp_points);
    iret = Util_TableGetIntArray (param_table_handle, N_interp_points,
                                  &source_map.front(), "source_map");
    have_source_map = not (iret == UTIL_ERROR_TABLE_NO_SUCH_KEY);
    if (not have_source_map) {
      // No explicit source map specified
      if (Carpet::map != -1) {
        // Interpolate from the current map
        source_map.assign (source_map.size(), Carpet::map);
      } else if (maps == 1) {
        // Interpolate from map 0 if this is the only one
        // (for backwards compatibility)
        source_map.assign (source_map.size(), 0);
      } else {
        CCTK_WARN (CCTK_WARN_ALERT, "No source map specified");
        return -1;
      }
    } else if (iret < 0) {
      CCTK_WARN (CCTK_WARN_ALERT, "internal error");
      return -1;
    } else if (iret != N_interp_points) {
      CCTK_WARN (CCTK_WARN_ALERT, "Source map array has wrong size");
      return -1;
    } else {
      iret = Util_TableGetIntArray (param_table_handle, source_map.size(),
                                    &source_map.front(), "source_map");
      assert (iret == (int)source_map.size());

#ifndef _NDEBUG
      // Check source map
#pragma omp parallel for
      for (int n = 0; n < (int)source_map.size(); ++n) {
        if (not (source_map[n] >= 0 and source_map[n] < maps)) {
          cout << "CI: n=" << n << " map=" << source_map[n] << endl;
        }
        assert (source_map[n] >= 0 and source_map[n] < maps);
      }
#endif
    }

    return 0;
  }



  // Find the component and integer index to which a grid point
  // belongs.  This uses a linear search over all components, which
  // does NOT scale with the number of components.
  static
  void
  find_location_linear (gh const * restrict const hh,
                        rvect const & restrict pos,
                        rvect const & restrict lower,
                        rvect const & restrict upper,
                        rvect const & restrict delta,
                        int const ml,
                        int const minrl, int const maxrl,
                        int & restrict rl,
                        int & restrict c)
  {
    // cout << "CarpetInterp: assign: m=" << m << " pos=" << pos << endl;

    assert (ml>=0 and ml<mglevels);
    assert (minrl>=0 and minrl<maxrl and maxrl<=reflevels);

    CCTK_REAL const rone  = 1.0;
    CCTK_REAL const rhalf = rone / 2;

    if (all (pos >= lower and pos <= upper)) {
      // The point is within the domain

      // Try finer levels first
      for (rl = maxrl-1; rl >= minrl; --rl) {

        ivect const fact =
          maxspacereflevelfact / spacereffacts.AT(rl) * ipow(mgfact, ml);
        ivect const ipos =
          ivect(floor((pos - lower) / (delta * rvect(fact)) + rhalf)) * fact;

        ivect const & stride = hh->baseextent(ml,rl).stride();
        assert (all (ipos % stride == 0));

        gh::cregs const & regs = hh->regions.AT(ml).AT(rl);

        // Search all components linearly
        for (c = 0; c < int(regs.size()); ++c) {
          region_t const & reg = regs.AT(c);
          if (reg.extent.contains(ipos)) {
            // We found the refinement level, component, and index to
            // which this grid point belongs
            return;
          }
        }
      }
    }

    // The point does not belong to any component.  This should happen
    // only rarely.
    rl = -1;
    c = -1;
  }



  // Find the component and integer index to which a grid point
  // belongs.  This uses a tree search over the superregions in the
  // grid struction, which should scale reasonably (O(n log n)) better
  // with the number of componets components.
  static
  void
  find_location_tree (gh const * restrict const hh,
                      rvect const & restrict pos,
                      rvect const & restrict lower,
                      rvect const & restrict upper,
                      rvect const & restrict delta,
                      int const ml,
                      int const minrl, int const maxrl,
                      int & restrict rl,
                      int & restrict c)
  {
    // cout << "CarpetInterp: assign: m=" << m << " pos=" << pos << endl;

    assert (ml>=0 and ml<mglevels);
    assert (minrl>=0 and minrl<maxrl and maxrl<=reflevels);

    CCTK_REAL const rone  = 1.0;
    CCTK_REAL const rhalf = rone / 2;

    if (all (pos >= lower and pos <= upper)) {
      // The point is within the domain

      // Try finer levels first
      for (rl = maxrl-1; rl >= minrl; --rl) {

        ivect const fact =
          maxspacereflevelfact / spacereffacts.AT(rl) * ipow(mgfact, ml);
        ivect const ipos =
          ivect(floor((pos - lower) / (delta * rvect(fact)) + rhalf)) * fact;

        ivect const & stride = hh->baseextent(ml,rl).stride();
        assert (all (ipos % stride == 0));

        gh::cregs const & regs = hh->superregions.AT(rl);

        // Search all superregions linearly.  Each superregion
        // corresponds to a "refined region", and the number of
        // superregions is thus presumably independent of the number
        // of processors.
        for (size_t r = 0; r < regs.size(); ++r) {
          region_t const & reg = regs.AT(r);
          if (reg.extent.contains(ipos)) {
            // We found the superregion to which this grid point
            // belongs

            // Search the superregion hierarchically
            pseudoregion_t const * const preg = reg.processors->search(ipos);
            assert (preg);

            // We now know the refinement level, component, and index
            // to which this grid point belongs
            c = preg->component;
            return;
          }
        }
      }
    }

    // The point does not belong to any component.  This should happen
    // only rarely.
    rl = -1;
    c = -1;
  }



  static void
  map_points (cGH const* const cctkGH,
              int const coord_system_handle,
              int const coord_group,
              int const ml,
              int const minrl,
              int const maxrl,
              int const maxncomps,
              int const N_dims,
              int const npoints,
              vector<CCTK_INT> & source_map,
              void const* const coords_list[],
              CCTK_INT * const procs,
              CCTK_INT * const rlev)
  {
    DECLARE_CCTK_PARAMETERS;

    static Timer * timer = NULL;
    if (not timer) timer = new Timer ("MapPoints::map_points");
    timer->start ();

    assert (npoints == 0 or coords_list);

    assert ((int)source_map.size() == npoints);


    // Find out about the coordinates: origin and delta for the Carpet
    // grid indices
    vector<rvect> lower (maps);
    vector<rvect> upper (maps);
    vector<rvect> delta (maps); // spacing on finest possible grid

    int const grouptype = CCTK_GroupTypeI (coord_group);
    assert(grouptype == CCTK_GF);
    for (int m = 0; m < Carpet::maps; ++ m) {
      jvect gsh;
      GetCoordRange (cctkGH, m, mglevel, dim,
                     & gsh[0],
                     & lower.AT(m)[0], & upper.AT(m)[0], & delta.AT(m)[0]);
      delta.AT(m) /= maxspacereflevelfact;
    }

    // Assign interpolation points to processors/components
#pragma omp parallel for
    for (int n = 0; n < npoints; ++n) {

      int & m = source_map.AT(n);
      int rl, c = -1;
      rvect pos;
      gh const * hh = NULL;

      if (m >= 0) {

        hh = arrdata.AT(coord_group).AT(m).hh;

        // Find the grid point closest to the interpolation point
        for (int d = 0; d < N_dims; ++d) {
          pos[d] = static_cast<CCTK_REAL const *>(coords_list[d])[n];
        }

        // Find the component that this grid point belongs to

        // Calculate rl, c, and proc
        if (not tree_search) {
          find_location_linear
            (hh, pos, lower.AT(m), upper.AT(m), delta.AT(m), ml, minrl, maxrl,
           rl, c);
        } else {
          find_location_tree
            (hh, pos, lower.AT(m), upper.AT(m), delta.AT(m), ml, minrl, maxrl,
             rl, c);
          if (check_tree_search) {
            int rl2, c2;
            find_location_linear
              (hh, pos, lower.AT(m), upper.AT(m), delta.AT(m), ml, minrl, maxrl,
               rl2, c2);
            if (rl2!=rl or c2!=c) {
#pragma omp critical
              CCTK_VWarn (CCTK_WARN_ABORT, __LINE__, __FILE__, CCTK_THORNSTRING,
                          "Inconsistent search result from find_location_tree for interpolation point #%d at [%g,%g,%g] of patch #%d is not on any component",
                          n, (double)pos[0], (double)pos[1], (double)pos[2], (int)m);
            }
          }
        }

      } // if m >= 0

      if (c == -1) {
        // The point could not be mapped onto any component

        // Warn only once, namely when mapping points onto processors.
        // (This routine is called twice; first to map points onto
        // processors, then to map points onto components.)
#pragma omp critical
        CCTK_VWarn (CCTK_WARN_PICKY, __LINE__, __FILE__, CCTK_THORNSTRING,
                    "Interpolation point #%d at [%g,%g,%g] of patch #%d is not on any component",
                    n, (double)pos[0], (double)pos[1], (double)pos[2], (int)m);

        // Map the point (arbitrarily) to the first component of the
        // coarsest grid
        // TODO: Handle these points explicitly later on
        rl = minrl;
        c = 0;

        // Find a patch which exists on this processor
        for (m=0; m<maps; ++m) {
          hh = arrdata.AT(coord_group).AT(m).hh;
          if (hh->components(rl) > c) break;
        }
        assert (m < maps);
      }

#ifndef _NDEBUG
      if (not (rl >= minrl and rl < maxrl) or
          not (c >= 0 and c < hh->components(rl)))
      {
        cout << "CI: m=" << m << " rl=" << rl << " c=" << c << " ext=" << hh->extent(ml,rl,c) << endl;
      }
      assert (rl >= minrl and rl < maxrl);
      assert (c >= 0 and c < hh->components(rl));
#endif

      int const proc = hh->processor(rl,c);
      procs[n] = proc;
      rlev[n] = rl;

    } // for n

    timer->stop (npoints);
  }
} 

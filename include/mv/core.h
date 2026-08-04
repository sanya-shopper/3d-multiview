#ifndef MV_CORE_H
#define MV_CORE_H

/* Paper: doc/multiview.tex, section Conventions and module layout: "System overview". */

/* multiview — portable C99 multiview-geometry library.
 *
 * Conventions used throughout:
 *  - all matrices are row-major double arrays; a "3x3" is double[9]
 *  - points are Euclidean unless a function says homogeneous
 *  - pixel coordinates: u right, v down, origin at the top-left corner
 *  - world-to-camera transform: Xc = R * X + t
 */

#define MV_VERSION "0.6.0"

enum { MV_OK = 0, MV_ERR = -1 };

#endif /* MV_CORE_H */

#ifndef model_h
#define model_h

#include "dtype.h"

/**
 * @brief Constructs a quantum link model for a triangular lattice.
 *
 * This function generates a quantum link model on a triangular lattice 
 * with given parameters. The lattice is composed of sublattice A and 
 * sublattice B, and the model includes various types of bonds:
 * temp_trans, reference_vonf, and gauss_law.
 *
 * @param[in] lx     The number of lattice sites in the x direction.
 * @param[in] ly     The number of lattice sites in the y direction.
 * @param[in] lambda The weight of the reference_vonf bonds.
 * @return A pointer to the constructed model.
 *
 * @note The generated model should be freed using the appropriate function 
 *  when no longer needed.
 */
model* quantum_link_triangular_lattice(int lx, int ly, double lambda);


#endif

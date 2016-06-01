/*
 * File created originally on Jul 7, 2005
 */
package nl.tudelft.bt.model.multigrid.boundary_conditions;

/**
 * Interface for the boundary conditions to be used in the multigrid method
 * 
 * @author jxavier
 */
public interface BoundaryConditions {
	/**
	 * Refreshes the boundary conditions during multigrid computatioin
	 * 
	 * @param u a 3D matrix to update boundary conditions in.
	 */
	public void refreshBoundaryConditions(float u[][][]);
	/**
	 * Implements the shape of the substratum carrier for the multigrid
	 * computation
	 * 
	 * @param i
	 * @param j
	 * @param k
	 * @return true if i,j,k is in the carrier, false otherwise
	 */
	public boolean isCarrier(int i, int j, int k);
}

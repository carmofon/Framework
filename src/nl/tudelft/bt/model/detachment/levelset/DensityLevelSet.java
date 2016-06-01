/*
 * Created on 29-jan-2004 by Joao Xavier (j.xavier@tnw.tudelft.nl)
 */
package nl.tudelft.bt.model.detachment.levelset;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import nl.tudelft.bt.model.ContinuousCoordinate;
import nl.tudelft.bt.model.Model;
import nl.tudelft.bt.model.detachment.GridElementFloat;
import nl.tudelft.bt.model.detachment.GridElementFloat.TValueComparator;
import nl.tudelft.bt.model.detachment.levelset.functions.DetachmentSpeedFunction;
import nl.tudelft.bt.model.exceptions.ModelRuntimeException;
import nl.tudelft.bt.model.multigrid.MultigridUtils;
import nl.tudelft.bt.model.particlebased.BiomassParticle;
import nl.tudelft.bt.model.particlebased.BiomassParticleContainer;
import nl.tudelft.bt.model.particlebased.BiomassParticleSet;
import nl.tudelft.bt.model.util.ExtraMath;

/**
 * Implements determination of level set for biomass erosion based on density
 * 
 * @author Joao Xavier (j.xavier@tnw.tudelft.nl)
 */
public class DensityLevelSet {
	private static final float INF = Float.POSITIVE_INFINITY;

	private int _n;

	private int _m;

	private int _l;

	private float _shovingGridSide;

	private DetachmentSpeedFunction _detachmentFunction;

	private float _timeStep;

	private float[][][] _levelset;

	private BiomassParticleSet[][][] _shovingGrid;

	private boolean[][][] _bulkLiquidElements;

	private boolean[][][] _carrierElements;

	private ArrayList _alive;

	private ArrayList _close;

	/**
	 * Create a new instance of DensityLevelSet and compute the actual level
	 * set, which can be recovered through method getLevelSetMatrix
	 * 
	 * @param carrierElements
	 *            TODO
	 */
	public DensityLevelSet(BiomassParticleSet[][][] shovingGrid,
			boolean[][][] bulkLiquidElements, boolean[][][] carrierElements,
			float shovingGridSide, DetachmentSpeedFunction detachmentFunction,
			float timeStep) {
		// assign
		_shovingGrid = shovingGrid;
		_bulkLiquidElements = bulkLiquidElements;
		_carrierElements = carrierElements;
		_shovingGridSide = shovingGridSide;
		_detachmentFunction = detachmentFunction;
		_timeStep = timeStep;
		//initialize containers
		_alive = new ArrayList();
		_close = new ArrayList();
		//get size of shoving grid
		_n = bulkLiquidElements.length;
		_m = bulkLiquidElements[0].length;
		_l = bulkLiquidElements[0][0].length;
		//create levelset holder
		_levelset = new float[_n][_m][_l];
		// go through the whole matrix
		//finding and initializing the border points
		for (int i = 0; i < _n; i++) {
			for (int j = 0; j < _m; j++) {
				for (int k = 0; k < _l; k++) {
					setBorderPointValue(i, j, k);
				}
			}
		}
		//main loop
		while (_close.size() > 0) {
			// Order the elemets in _close in respect to the T vlaue
			Collections.sort(_close, new GridElementFloat.TValueComparator());
			// get the first element (trial)
			GridElementFloat trial = (GridElementFloat) _close.get(0);
			// remove trial from _close and put in _alive
			_close.remove(trial);
			_alive.add(trial);
			// get all neighbors of trial that do not belog to _close
			int i = trial.i;
			int j = trial.j;
			int k = trial.k;
			// bottom neighbor
			if (i > 0)
				addToCloseAndUpdate(i - 1, j, k);
			// y-side neighbors:
			addToCloseAndUpdate(i, (j < _m - 1) ? (j + 1) : 0, k);
			addToCloseAndUpdate(i, (j > 0) ? (j - 1) : (_m - 1), k);
			// z-side neighbors:
			if (_l > 1) {
				addToCloseAndUpdate(i, j, (k < _l - 1) ? k + 1 : 0);
				addToCloseAndUpdate(i, j, (k > 0) ? (k - 1) : (_l - 1));
			}
			// top neighbor:
			if (i < _n - 1)
				addToCloseAndUpdate(i + 1, j, k);
		}
	}

	/**
	 * In case the element is not empty and does not belong to _close, add it
	 * and compute the t value
	 * 
	 * @param i
	 * @param j
	 * @param k
	 */
	private void addToCloseAndUpdate(int i, int j, int k) {
		// check if this is a point to compute
		if (!_bulkLiquidElements[i][j][k] & (!isEmptyCarrier(i, j, k))) {
			GridElementFloat element = new GridElementFloat(_levelset, i, j, k);
			if (!_close.contains(element) & !_alive.contains(element)) {
				// compute the T value
				float t = computeTValue(i, j, k);
				// make sure t is a valid number
				if (Float.isNaN(t) | (t < 0)) {
					throw new ModelRuntimeException(
							"Illegal T value computed in levelset " + t);
				}
				_levelset[i][j][k] = t;
				// mark all particles for detachment if its the case
				if (t < _timeStep)
					markParticlesForErosion(i, j, k);
				//add element to _clos
				_close.add(element);
			}
		}
	}

	/**
	 * Determine the level set for detachment based on biomass density
	 * 
	 * @param shovingGrid
	 *            the shoving grid
	 * @param detachmentRate
	 *            the detachment rate
	 * @return the levelset
	 */
	private void setBorderPointValue(int i, int j, int k) {
		if (_bulkLiquidElements[i][j][k])
			_levelset[i][j][k] = 0; // bulk liquid
		else {
			_levelset[i][j][k] = INF;
			if (!isEmptyCarrier(i, j, k)) {
				// only determine the point if this is not carrier
				int nn = numberOfBulkLiquidNeighbors(i, j, k);
				if (nn > 0) { //Border points
					_close.add(new GridElementFloat(_levelset, i, j, k));
					// The time necessary for full erosion of a volume ement at
					// the
					// border based on its mass and exposed area
					float t = _shovingGridSide
							/ (getLocalDetachmentSpeed(i, j, k) * nn);
					// make sure t is a valid number
					if (Float.isNaN(t) | (t < 0))
						throw new ModelRuntimeException(t);
					// and update the levelset
					_levelset[i][j][k] = t;
					// mark all particles for detachment if its the case
					if (t < _timeStep)
						markParticlesForErosion(i, j, k);
				}
			}
		}
	}

	/**
	 * Set all particles in the grid elemente to detach
	 * 
	 * @param i
	 * @param j
	 * @param k
	 */
	private void markParticlesForErosion(int i, int j, int k) {
		BiomassParticleSet ps = _shovingGrid[i][j][k];
		if (ps != null) {
			ps.resetIndex();
			while (ps.hasNext()) {
				// add particle to detached biomass and mark it for detachmnent
				BiomassParticle p = ps.next();
				p.setToDetachByErosion();
			}
		}
	}

	/**
	 * The number of neighbors that containe no biomass
	 * 
	 * @param l
	 * @param m
	 * @param n
	 * @return the number of empty neighbors
	 */
	private int numberOfBulkLiquidNeighbors(int i, int j, int k) {
		int nn = 0;
		// x neighbors
		// bottom neighbor
		if (i > 0)
			if (isBulkLiquid(i - 1, j, k))
				nn++;
		// top neighbor:
		if (i < _n - 1)
			if (isBulkLiquid(i + 1, j, k))
				nn++;
		// y-side neighbors:
		if (isBulkLiquid(i, j < _m - 1 ? j + 1 : 0, k))
			nn++;
		if (isBulkLiquid(i, j > 0 ? j - 1 : _m - 1, k))
			nn++;
		// z-side neighbors:
		if (_l > 1) {
			if (isBulkLiquid(i, j, k < _l - 1 ? k + 1 : 0))
				nn++;
			if (isBulkLiquid(i, j, k > 0 ? k - 1 : _l - 1))
				nn++;
		}
		return nn;
	}

	/**
	 * Return true if the elemnt belongs to bulk liquid
	 * 
	 * @param n
	 * @param l
	 * @param m
	 * @return true if the elemnt belongs to bulk liquid
	 */
	private boolean isBulkLiquid(int i, int j, int k) {
		return _bulkLiquidElements[i][j][k];
	}

	/**
	 * Check if entry i, j, k in the shoving grid is empty carrier (carrier
	 * elements at the border might also have biomass)
	 * 
	 * @param i
	 * @param j
	 * @param k
	 * @return true if entry i, j, k has no particles
	 */
	private boolean isEmptyCarrier(int i, int j, int k) {
		BiomassParticleSet ps = _shovingGrid[i][j][k];
		if (ps != null) {
			// return false if there are particles in the system
			if (ps.getNumberOfParticles() > 0)
				return false;
		}
		// if this has no biomass, return true only if there is
		// no liquid and some carrier
		return _carrierElements[i][j][k];
	}

	/**
	 * Get the new T value for level set as based on the values of the neighbors
	 * and marks the particle for detachment if t is less than the time step
	 * 
	 * @param i
	 * @param j
	 * @param k
	 * @return next levelset value
	 */
	private float computeTValue(int i, int j, int k) {
		//X (does not compute the bottom neighbor for lower points)
		float tXminus = (i > 0 ? _levelset[i - 1][j][k] : INF);
		float tXplus = (i < _n ? _levelset[i + 1][j][k] : INF);
		float tX = Math.min(tXminus, tXplus);
		//Y (implements periodicity)
		float tYminus = _levelset[i][j > 0 ? j - 1 : _m - 1][k];
		float tYplus = _levelset[i][j < _m - 1 ? j + 1 : 0][k];
		float tY = Math.min(tYminus, tYplus);
		//Z (implements periodicity)
		// if 2D, z values will be the same and equal to present element value
		float tZminus = _levelset[i][j][k > 0 ? k - 1 : _l - 1];
		float tZplus = _levelset[i][j][k < _l - 1 ? k + 1 : 0];
		float tZ = Math.min(tZminus, tZplus);
		// compute the solution for all possible combinations and choose the
		// maximum value among the valid solutions
		float validSolution = 0;
		float approximateSolution = 0;
		// each difference value will be tried once in any combination
		for (float fi = tX;; fi = INF) {
			for (float fj = tY;; fj = INF) {
				for (float fk = tZ;; fk = INF) {
					// skip iteration where all are INF
					if (fi == INF && fj == INF && fk == INF)
						break;
					// get the roots
					float tmp = computeRoots(fi, fj, fk,
							getLocalDetachmentSpeed(i, j, k));
					// if tmp is a number, compute maximum for approximate
					// solution, else keep aproximate solution
					approximateSolution = (Float.isNaN(tmp)
							? approximateSolution
							: Math.max(tmp, approximateSolution));
					// check if solution is valid
					if (!Float.isNaN(tmp)) {
						approximateSolution = Math
								.max(tmp, approximateSolution);
						if (solutionNotValid(tmp, fi, tX)
								| solutionNotValid(tmp, fj, tY)
								| solutionNotValid(tmp, fk, tZ)) {
							// solution not valid, do nothing
						} else {
							//if flow reaches this point, solution is valid.
							// confront with previous maximum solution
							validSolution = Math.max(tmp, validSolution);
						}
					}
					// break
					if (fk == INF)
						break;
				}
				// break
				if (fj == INF)
					break;
			}
			// break
			if (fi == INF)
				break;
		}
		// validity check may return invalid in special cases due to precision
		// of float computations. In these cases, the approximate solution is
		// used
		return (validSolution == 0 ? approximateSolution : validSolution);
	}

	/**
	 * compute the maximum value of roots
	 * 
	 * @param tX
	 * @param tY
	 * @param tZ
	 * @param detachmentRate
	 * @return
	 */
	private float computeRoots(float tX, float tY, float tZ,
			float detachmentRate) {
		// parameters for solving quadratic equation
		float a = (tX != INF ? 1.0f : 0) + (tY != INF ? 1.0f : 0)
				+ (tZ != INF ? 1.0f : 0);
		float b = -2.0f
				* ((tX != INF ? tX : 0) + (tY != INF ? tY : 0) + (tZ != INF
						? tZ
						: 0));
		float c = (tX != INF ? tX * tX : 0) + (tY != INF ? tY * tY : 0)
				+ (tZ != INF ? tZ * tZ : 0)
				- ExtraMath.sq(_shovingGridSide / detachmentRate);
		// get the 2 solutions
		float aux = (float) Math.sqrt(b * b - 4.0f * a * c);
		//Positive solution is always the valid one
		return (-b + aux) / (2.0f * a);
	}

	/**
	 * Check the validity of a solution
	 * 
	 * @param s
	 *            solution value
	 * @param f
	 *            the present f (t or INF)
	 * @param t
	 *            the present t (min(tplus, tminus)
	 * @return true if solution is not valid
	 */
	private boolean solutionNotValid(float s, float f, float t) {
		//check validity criteria
		return ((f != INF) ? (s < t) : (s > t));
	}

	/**
	 * Set the bulk liquid elements matrix
	 * 
	 * @param liquidElements
	 *            The _bulkLiquidElements to set.
	 */
	public void setBulkLiquidElements(boolean[][][] liquidElements) {
		_bulkLiquidElements = liquidElements;
	}

	/**
	 * Performs erosion of particles that reamin at the border
	 * 
	 * @throws ModelException
	 *             if fraction to erode has illegal value
	 */
	public void erodeBorder() {
		//find the border points
		for (int i = 0; i < _n; i++) {
			for (int j = 0; j < _m; j++) {
				for (int k = 0; k < _l; k++) {
					int nn = numberOfBulkLiquidNeighbors(i, j, k);
					if (nn > 0) { //Border points
						// The time necessary for full erosion of a volume
						// ement at the
						// border based on its mass and exposed area
						float t = _levelset[i][j][k];
						BiomassParticleSet ps = _shovingGrid[i][j][k];
						if (ps != null) {
							ps.resetIndex();
							while (ps.hasNext()) {
								// erode particle and add detached biomass to
								// container
								Model.model().detachedBiomassContainer()
										.addToErodedBiomassSeries(
												ps.next().erode(_timeStep / t));
							}
						}
					}
				}
			}
		}
	}

	/**
	 * @param i
	 * @param j
	 * @param k
	 * @return the detachment speed at location [um.h^-1]
	 */
	private float getLocalDetachmentSpeed(int i, int j, int k) {
		ContinuousCoordinate c = ((BiomassParticleContainer) Model.model()
				.getBiomassContainer()).getGridElementCenter(i, j, k);
		return _detachmentFunction.getValue(c);
	}

	/**
	 * Write the current level set matrix to file
	 * 
	 * @param s
	 *            the filename
	 * @throws if
	 *             writing is not possible
	 */
	public void writeToFile(String s) throws IOException {
		File f = new File(s + ".txt");
		java.io.FileWriter fr = new java.io.FileWriter(f);
		fr.write(MultigridUtils.matrixToString(_levelset));
		fr.close();
	}
}
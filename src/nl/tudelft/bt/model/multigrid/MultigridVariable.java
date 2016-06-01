package nl.tudelft.bt.model.multigrid;

import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.Serializable;
import java.util.Collection;
import java.util.Iterator;

import nl.tudelft.bt.model.*;
import nl.tudelft.bt.model.exceptions.*;
import nl.tudelft.bt.model.multigrid.boundary_conditions.BoundaryConditions;
import nl.tudelft.bt.model.multigrid.boundary_layers.BoundaryLayer;

/**
 * Implements a 3D state variable at multigrid resolutions. Common base class
 * for BacteriaSpecies and ChemicalSpecies and implementation for relative
 * diffusion. Boundary conditions are implemented by a mechanism of padding
 * elements. When values of variable are changed, the boundary conditions are
 * updated
 * 
 * @see nl.tudelft.bt.model.multigrid.ParticulateSpecies
 * @see nl.tudelft.bt.model.multigrid.SoluteSpecies
 * @author Joao Xavier (j.xavier@tnw.tudelft.nl)
 */
public class MultigridVariable implements Serializable {
	// boundarylayer threshold
	protected static final float BLTHRESH = 0.5f;

	// index of coarsestgrid
	protected static final int COARSEST = 0;

	// Default values. May be reset using setSteps()
	protected static final int VCYCLES = 10;

	protected static final int NSOLVE = 1000; // coarset grid solution (50

	// works)

	protected static int nPreSteps = 50;

	// pre-smoothing steps (5 works on simple systems, 50 works on double
	// substrate limitation)
	protected static int nPosSteps = 500;

	// post-smoothing steps (50 works on simple systems, 500 works on double
	// substrate limitation)
	protected static final float ALPHA = 0.33f; // post-smoothing steps

	private static final float PRECISION = 1e-6f; // required precision

	protected String _name;

	// holds the multi-grid values for a variable
	// set access control to public when boundary layers where placed on
	// a different package
	public float[][][][] _mg;

	// sizes of finest grid
	protected static int _l;

	protected static int _m;

	protected static int _n; // vertical size

	// multigrid order
	protected static int _order;

	// the volume of an element of the finest grid [um^3]
	protected static float _voxelVolume;

	protected static float _voxelSide;

	protected static float _referenceSystemSide;

	// current grid order
	protected static int _g;

	// current indexes
	protected static int _i;

	protected static int _j;

	protected static int _k;
	
	// store the present location
	private static ContinuousCoordinate _presentLocation;

	// grid is set lock
	private static boolean _gridIsSet = false;

	// a static 3D array for neighborhood computations in relax
	private static final float[][][] _diff = new float[3][3][3];

	// a static array for rate and rate derivative computations in relax
	private static final float[] _rDr = new float[2];

	// an auxiliry variable
	private static final DiscreteCoordinate _auxDiscreteCoordinate = new DiscreteCoordinate();

	// an auxiliry variable
	private static final ContinuousCoordinate _auxContinuousCoordinate = new ContinuousCoordinate();

	// the boundary conditions
	protected static BoundaryConditions _boundaryConditions;

	// the reference to model (static)
	protected static Model _model = Model.model();

	/**
	 * Allocates memory for a multi-grid concentration matrix with padding
	 * entries for boundaries. Only valid for 3D cubic and 2D square systems.
	 * 
	 * @throws MultigridSystemNotSetException
	 */
	public MultigridVariable() throws MultigridSystemNotSetException {
		if (!_gridIsSet)
			throw new MultigridSystemNotSetException();
		_mg = new float[_order][][][];
		for (int i = 0; i < _order; i++) {
			int n = MultigridUtils.coarserSize(_n, _order - i - 1);
			int m = MultigridUtils.coarserSize(_m, _order - i - 1);
			int l = MultigridUtils.coarserSize(_l, _order - i - 1);
			// with padding for boundary conditions
			_mg[i] = new float[n + 2][m + 2][l + 2];
		}
	}

	/**
	 * Set the discrete grid properties. The grid size must be of the for n =
	 * 2^i + 1, where i is an integer. This implementation is for 3D cubic and
	 * 2D square systems.
	 * 
	 * @param n
	 *            grid size
	 * @throws InvalidValueException
	 *             in case n is not of the form n = 2^i + 1
	 */
	static public void setGrid(int n, int m, int l)
			throws InvalidValueException {
		// compute orders and check validity of grid
		int on = MultigridUtils.order(n);
		int om = MultigridUtils.order(m);
		// determine the order of the system
		_order = Math.min(on, om);
		if (l != 1) {
			// 3D
			int ol = MultigridUtils.order(l);
			_order = Math.min(_order, ol);
		}
		// compute the volume of one finest grid element
		_referenceSystemSide = _model.referenceSystemSide;
		_voxelSide = _referenceSystemSide / m;
		_voxelVolume = (float) Math.pow(_voxelSide, 3);
		// finest grid size
		_n = n;
		_m = m;
		_l = l;
		// open the lock, so new variables can be created
		_gridIsSet = true;
	}

	/**
	 * Convert continuous coordinate into indexes of the discrete grid
	 * 
	 * @param c
	 *            coordinate to converted
	 * @return contains corresponding indexes in the discrete grid
	 */
	public static DiscreteCoordinate snapToDiscrete(ContinuousCoordinate c) {
		// snap coordinate to grid indexes
		_auxDiscreteCoordinate.i = snapToDiscreteI(c.x);
		_auxDiscreteCoordinate.j = snapToDiscreteJ(c.y);
		_auxDiscreteCoordinate.k = snapToDiscreteJ(c.z);
		return _auxDiscreteCoordinate;
	}

	public static int snapToDiscreteI(float x) {
		// NOTE + 1 acounts for padding
		return (int) ((x / _model.systemSize.x) * (float) _n) + 1;
	}

	public static int snapToDiscreteJ(float y) {
		// NOTE + 1 acounts for padding
		return (int) ((y / _model.systemSize.y) * (float) _m) + 1;
	}

	public static int snapToDiscreteK(float z) {
		if (z > _model.systemSize.z)
			// this code is specific for the 2D case
			return 1;
		else
			// this code is specific for the 3D case
			return (int) ((z / _model.systemSize.z) * (float) _l) + 1;
		// NOTE + 1 acounts for padding
	}

	/**
	 * Convert a discrete coordinate (indexes of matrix in finest grid) into
	 * Cartesian coordinates
	 * 
	 * @param d
	 * @return the Cartesian location
	 */
	public static ContinuousCoordinate convertDiscreteToContinuous(
			DiscreteCoordinate d) {
		// add 0.5 to the coordinated so that the coordinate returned
		// is the center of the grid node
		_auxContinuousCoordinate.x = convertDiscreteToContinuousX(d.i);
		_auxContinuousCoordinate.y = convertDiscreteToContinuousY(d.j);
		_auxContinuousCoordinate.z = convertDiscreteToContinuousZ(d.k);
		return _auxContinuousCoordinate;
	}

	public static float convertDiscreteToContinuousX(int i) {
		return (((float) i) - 0.5f) * _model.systemSize.x / ((float) _n);
	}

	public static float convertDiscreteToContinuousY(int j) {
		return (((float) j) - 0.5f) * _model.systemSize.y / ((float) _m);
	}

	public static float convertDiscreteToContinuousZ(int k) {
		return (((float) k) - 0.5f) * _model.systemSize.z / ((float) _l);
	}

	/**
	 * Set the value of the present location for retrieval of values
	 * 
	 * @param c
	 */
	public static void setCurrentLocation(ContinuousCoordinate c) {
		_presentLocation = c;
		DiscreteCoordinate d = snapToDiscrete(c);
		_i = d.i;
		_j = d.j;
		_k = d.k;
	}

	
	/**
	 * Get the value of the present location
	 * 
	 * @return the present location set by multigrid calculation
	 */
	public static ContinuousCoordinate getCurrentLocation() {
		return _presentLocation;
	}

	
	/**
	 * Get the current value
	 * 
	 * @return value for current multigrid entry
	 */
	public float getValue() {
		return _mg[_g][_i][_j][_k];
	}

	/**
	 * Gets the value of concentration of this chemical at location c
	 * 
	 * @param c
	 *            location to get value
	 * @return concentration at c
	 */
	public float getValueAt(ContinuousCoordinate c) {
		DiscreteCoordinate dc = snapToDiscrete(c);
		return _mg[_order - 1][dc.i][dc.j][dc.k];
	}

	/**
	 * Gets the value of concentration of this chemical at location i, j, k
	 * 
	 * @param i
	 * @param j
	 * @param k
	 * @return value at location i, j, k in finest grid
	 */
	public float getValueAt(int i, int j, int k) {
		return _mg[_order - 1][i][j][k];
	}

	/**
	 * Get the maximum value for variable in finest grid
	 * 
	 * @return maximum concentration of solute in system
	 */
	public float getMaximumValue() {
		return MultigridUtils.max(_mg[_order - 1]);
	}

	/**
	 * Get the minimum value for variable in finest grid
	 * 
	 * @return minimum concentration of solute in system
	 */
	public float getMinimumValue() {
		return MultigridUtils.min(_mg[_order - 1]);
	}

	/**
	 * Increment value at finest grid at position c by the value f
	 * 
	 * @param c
	 *            position to increment value
	 * @param f
	 *            value to incrment
	 */
	void incrementValueAt(ContinuousCoordinate c, float f) {
		DiscreteCoordinate dc = snapToDiscrete(c);
		_mg[_order - 1][dc.i][dc.j][dc.k] += f;
	}

	/**
	 * Set all entries of the finest grid to a given value
	 * 
	 * @param val
	 *            value to set entries to
	 */
	protected void setValue(float val) {
		MultigridUtils.setValues(_mg[_order - 1], val);
	}

	/**
	 * Update the multigrid data by restricting the finner grid progrssivelly
	 * down to the coarsest grid.
	 */
	protected void updateMultigridCopies() {
		for (int i = _order - 1; i > COARSEST; i--) {
			MultigridUtils.restrict(_mg[i], _mg[i - 1], _boundaryConditions);
		}
	}

	/**
	 * Solve the reaction and diffusion equations to pseudo- steady state. As a
	 * final result, the concentrations of chemical species are updated. The
	 * algorithm starts by udating the matrices of discrete data for the
	 * bacteria species, after which it solves the rection diffusion equations
	 * for all chemical sepcies using multigrid.
	 * 
	 * @param chem
	 * @param bac
	 * @param bl
	 * @throws MultigridSystemNotSetException
	 * 
	 */
	public static void solve(SoluteSpecies[] chem, ParticulateSpecies[] bac,
			BoundaryLayer bl) throws MultigridSystemNotSetException {
		// two temporary multigrid variables are needed for the
		// computation, are initialized with 0 allover
		MultigridVariable itemp = new MultigridVariable();
		MultigridVariable itau = new MultigridVariable();
		// create a relative diffusion data, so far it is 1 everywhere
		RelativeDiffusion relDiff = new RelativeDiffusion();
		bl.setBoundaryLayer(bac, _boundaryConditions);
		relDiff.computeValues(bac, bl, _boundaryConditions);
		// create multigrid copies of biomass and relative
		// diffusion
		for (int i = 0; i < bac.length; i++) {
			bac[i].updateMultigridCopies();
		}
		relDiff.updateMultigridCopies();
		bl.updateMultigridCopies();
		// Initialize concentration of all chemicals to value of
		// bulk concentration
		for (int i = 0; i < chem.length; i++) {
			chem[i].resetMultigridCopies();
		}
		// solve chemical concentrations on coarsest grid
		solveCoarsest(chem, relDiff, bl);
		// nested iteration loop
		for (int outer = 1; outer < _order; outer++) {
			_g = outer;
			for (int i = 0; i < chem.length; i++) {
				MultigridUtils.interpolateBoundaryLayer(chem[i]._mg[_g],
						chem[i]._mg[_g - 1], bl._mg[_g], _boundaryConditions);
				// set each chemical's r.h.s. to 0
				MultigridUtils.setValues(chem[i].rhs._mg[_g], 0.0f);
			}
			// V-cycle loop
			for (int v = 0; v < VCYCLES; v++) {
				// downward stroke of V
				while (_g > 0) {
					// pre-smoothing
					for (int j = 0; j < nPreSteps; j++) {
						relax(chem, relDiff, bl);
					}
					for (int j = 0; j < chem.length; j++) {
						// restrict uh to uH
						MultigridUtils.restrictBoundaryLayer(chem[j]._mg[_g],
								chem[j]._mg[_g - 1], bl._mg[_g - 1],
								_boundaryConditions);
						//
						lop(itemp, chem[j], relDiff, bl);
						//
						MultigridUtils.restrictBoundaryLayer(itemp._mg[_g],
								itemp._mg[_g - 1], bl._mg[_g - 1],
								_boundaryConditions);
						// reduce grid value _g temporarily
						_g--;
						lop(itau, chem[j], relDiff, bl);
						MultigridUtils.subtractTo(itau._mg[_g], itemp._mg[_g]);
						// sum tau to rhs of _g - 1
						MultigridUtils.restrictBoundaryLayer(
								chem[j].rhs._mg[_g + 1], chem[j].rhs._mg[_g],
								bl._mg[_g], _boundaryConditions);
						MultigridUtils.addTo(chem[j].rhs._mg[_g], itau._mg[_g]);
						// compute the truncation error for this V-cycle
						// for all chemicals
						if (_g + 1 == outer)
							chem[j].truncationError = ALPHA
									* MultigridUtils.computeNorm(itau._mg[_g]);
						// put _g value back for remaining solutes
						_g++;
					}
					// reduce grid value _g for good
					_g--;
				}
				// bottom of V
				solveCoarsest(chem, relDiff, bl);
				// upward stroke of V
				while (_g < outer) {
					_g++;
					for (int j = 0; j < chem.length; j++) {
						MultigridUtils.restrictBoundaryLayer(chem[j]._mg[_g],
								itemp._mg[_g - 1], bl._mg[_g - 1],
								_boundaryConditions);
						MultigridUtils.subtractTo(chem[j]._mg[_g - 1],
								itemp._mg[_g - 1]);
						MultigridUtils.interpolateBoundaryLayer(itau._mg[_g],
								chem[j]._mg[_g - 1], bl._mg[_g],
								_boundaryConditions);
						MultigridUtils.addTo(chem[j]._mg[_g], itau._mg[_g]);
					}
					truncatePresentSoluteValuesToZero(chem, bl);
					// post-smoothing
					for (int j = 0; j < nPosSteps; j++) {
						relax(chem, relDiff, bl);
					}
				}
				// break the V-cycles if remaining error is dominated
				// by local truncation error (see p. 884 of Numerical Recipes)
				boolean breakVCycle = true;
				for (int i = 0; i < chem.length; i++) {
					// compute the residue for this solute species
					lop(itemp, chem[i], relDiff, bl);
					MultigridUtils.subtractTo(itemp._mg[_g],
							chem[i].rhs._mg[_g]);
					float res = MultigridUtils.computeNorm(itemp._mg[_g]);
					// confirm that criterium is met for each solute
					if (v > 0) {
						System.out.println("grid " + _g + "; v " + v + "; "
								+ chem[i]._name + " res " + res + "; truncerr "
								+ chem[i].truncationError);
					}
					// confirm that criterium is met for each solute
					if (res > chem[i].truncationError) {
						breakVCycle = false;
						break;
					}
				}
				if (breakVCycle)
					break;
			}
		}
	}

	/**
	 * Find solution for the coarsest grid. Sets the current grid to coarsest,
	 * solutes values to 0 and relaxes NSOLVE times.
	 */
	private static void solveCoarsest(SoluteSpecies[] c, RelativeDiffusion d,
			BoundaryLayer bl) {
		_g = COARSEST;
		float err = 0;
		// reset coarsest grid to bulk concentration
		for (int i = 0; i < c.length; i++) {
			c[i].setValueCoarsestToBulk();
		}
		// relax NSOLVE times
		for (int j = 0; j < NSOLVE; j++) {
			relax(c, d, bl);
		}
	}

	/**
	 * Make sure all solute concentration values for the present are not lower
	 * than 0
	 * 
	 * @param c
	 *            array of the solute species in the system
	 * @param bl
	 *            the boundary layer
	 */
	private static void truncatePresentSoluteValuesToZero(SoluteSpecies[] c,
			BoundaryLayer bl) {
		int n = c[0]._mg[_g].length - 2;
		int m = c[0]._mg[_g][0].length - 2;
		int l = c[0]._mg[_g][0][0].length - 2;
		float v;
		for (_i = 1; _i <= n; _i++) {
			for (_j = 1; _j <= m; _j++) {
				for (_k = 1; _k <= l; _k++) {
					for (int chem = 0; chem < c.length; chem++) {
						if (bl._mg[_g][_i][_j][_k] < BLTHRESH) {
							v = c[chem]._mg[_g][_i][_j][_k];
							c[chem]._mg[_g][_i][_j][_k] = (v < 0 ? 0 : v);
						}
					}
				}
			}
		}
	}

	/**
	 * Perform relaxation for concentration of cehmical species at the current
	 * grid order
	 * 
	 * @param c
	 * @param d
	 * @param bl
	 */
	private static void relax(SoluteSpecies[] c, RelativeDiffusion d,
			BoundaryLayer bl) {
		float r, dr;
		int n = d._mg[_g].length - 2;
		int m = d._mg[_g][0].length - 2;
		int l = d._mg[_g][0][0].length - 2;
		float h = _referenceSystemSide / ((float) n - 1);
		float h2i = 0.5f / (h * h);
		// red-black relaxation
		// iterate through system
		// isw, jsw and ksw alternate between values 1 and 2
		int isw = 1;
		int jsw, ksw;
		for (int pass = 1; pass <= 2; pass++, isw = 3 - isw) {
			jsw = isw;
			for (_i = 1; _i <= n; _i++, jsw = 3 - jsw) {
				ksw = jsw;
				for (_j = 1; _j <= m; _j++, ksw = 3 - ksw) {
					for (_k = ksw; _k <= l; _k += 2) {
						for (int chem = 0; chem < c.length; chem++) {
							if (bl._mg[_g][_i][_j][_k] < BLTHRESH) {
								// Case: Inside boundary layer
								// Equations must be solved here
								float[][][] u = c[chem]._mg[_g];
								//
								c[chem]
										.updateValuesForRateAndRateDerivative(_rDr);
								r = _rDr[0];
								dr = _rDr[1];
								// compute diffusivity values
								// and that of surrounding neighbors
								float dc = c[chem].getDiffusivity();
								_diff[0][1][1] = dc * d._mg[_g][_i - 1][_j][_k];
								_diff[2][1][1] = dc * d._mg[_g][_i + 1][_j][_k];
								_diff[1][0][1] = dc * d._mg[_g][_i][_j - 1][_k];
								_diff[1][2][1] = dc * d._mg[_g][_i][_j + 1][_k];
								_diff[1][1][0] = dc * d._mg[_g][_i][_j][_k - 1];
								_diff[1][1][2] = dc * d._mg[_g][_i][_j][_k + 1];
								_diff[1][1][1] = dc * d._mg[_g][_i][_j][_k];
								// compute L operator
								float lop = ((_diff[2][1][1] + _diff[1][1][1])
										* (u[_i + 1][_j][_k] - u[_i][_j][_k])
										+ (_diff[0][1][1] + _diff[1][1][1])
										* (u[_i - 1][_j][_k] - u[_i][_j][_k])
										+ (_diff[1][2][1] + _diff[1][1][1])
										* (u[_i][_j + 1][_k] - u[_i][_j][_k])
										+ (_diff[1][0][1] + _diff[1][1][1])
										* (u[_i][_j - 1][_k] - u[_i][_j][_k])
										+ (_diff[1][1][2] + _diff[1][1][1])
										* (u[_i][_j][_k + 1] - u[_i][_j][_k]) + (_diff[1][1][0] + _diff[1][1][1])
										* (u[_i][_j][_k - 1] - u[_i][_j][_k]))
										* h2i + r;
								// compute derivative of L operator
								float dlop = -h2i
										* (6.0f * _diff[1][1][1]
												+ _diff[2][1][1]
												+ _diff[0][1][1]
												+ _diff[1][2][1]
												+ _diff[1][0][1]
												+ _diff[1][1][2] + _diff[1][1][0])
										+ dr;
								// compute residual
								float res = (lop - c[chem].rhs._mg[_g][_i][_j][_k])
										/ dlop;
								// update concentration (test for NaN)
								if (res != res) {
									System.out.println("---------------------");
									System.out.println("Multigrid problem:");
									System.out.println("_g = " + _g);
									System.out.println("_i = " + _i);
									System.out.println("_j = " + _j);
									System.out.println("_k = " + _k);
									System.out.println("r = " + r);
									System.out.println("dr = " + dr);
									System.out.println("concentration rhs = "
											+ c[chem].rhs._mg[_g][_i][_j][_k]);
									System.out.println("dlop = " + dlop);
									System.out
											.println("rhs/dlop= "
													+ (c[chem].rhs._mg[_g][_i][_j][_k] / dlop));
									System.out.println("solute species:");
									for (int s = 0; s < c.length; s++) {
										System.out.println(c[s]._name
												+ " (local) = "
												+ c[s]._mg[_g][_i][_j][_k]);
										System.out.println(c[s]._name
												+ " (bulk) = "
												+ c[s].getBulkConcentration());
									}
									// particulates:
									System.out.println("particulate species:");
									Collection ps = Model.model()
											.getParticulateSpecies();
									for (Iterator iter = ps.iterator(); iter
											.hasNext();) {
										ParticulateSpecies p = (ParticulateSpecies) iter
												.next();
										System.out.println(p._name + " = "
												+ p._mg[_g][_i][_j][_k]);
									}
									System.out.println("---------------------");
									String str = "NaN generated in multigrid solver "
											+ "while computing rate for "
											+ c[chem]._name;
									throw new ModelRuntimeException(str);
								}
								u[_i][_j][_k] -= res;
								// if negative concentrations, put 0 value
								u[_i][_j][_k] = (u[_i][_j][_k] < 0 ? 0
										: u[_i][_j][_k]);
							}
						}
					}
				}
			}
			// refresh the padding elements to enforce
			// boundary conditions for all solutes
			for (int i = 0; i < c.length; i++)
				_boundaryConditions.refreshBoundaryConditions(c[i]._mg[_g]);
		}
	}

	/**
	 * Compute the L-operator
	 * 
	 * @param res
	 * @param c
	 * @param d
	 * @param bl
	 *            boundary layer definition
	 */
	private static void lop(MultigridVariable res, SoluteSpecies c,
			RelativeDiffusion d, BoundaryLayer bl) {
		int n = d._mg[_g].length - 2;
		int m = d._mg[_g][0].length - 2;
		int l_ = d._mg[_g][0][0].length - 2;
		float h = _referenceSystemSide / ((float) n - 1);
		float h2i = 0.5f / (h * h);
		float lop; // temporary variable for L-operator
		// iterate through system
		for (_k = 1; _k <= l_; _k++) {
			for (_j = 1; _j <= m; _j++) {
				for (_i = 1; _i <= n; _i++)
					// compute lop only inside boundary layer
					if (bl._mg[_g][_i][_j][_k] < BLTHRESH) {
						// for simplification and easier access to
						// the current solute data:
						float[][][] u = c._mg[_g];
						// current rate for this solute
						float r = c.getRate();
						// compute diffusivity values
						// and that of surrounding neighbors
						float dc = c.getDiffusivity();
						_diff[0][1][1] = dc * d._mg[_g][_i - 1][_j][_k];
						_diff[2][1][1] = dc * d._mg[_g][_i + 1][_j][_k];
						_diff[1][0][1] = dc * d._mg[_g][_i][_j - 1][_k];
						_diff[1][2][1] = dc * d._mg[_g][_i][_j + 1][_k];
						_diff[1][1][0] = dc * d._mg[_g][_i][_j][_k - 1];
						_diff[1][1][2] = dc * d._mg[_g][_i][_j][_k + 1];
						_diff[1][1][1] = dc * d._mg[_g][_i][_j][_k];
						// compute L operator
						lop = ((_diff[2][1][1] + _diff[1][1][1])
								* (u[_i + 1][_j][_k] - u[_i][_j][_k])
								+ (_diff[0][1][1] + _diff[1][1][1])
								* (u[_i - 1][_j][_k] - u[_i][_j][_k])
								+ (_diff[1][2][1] + _diff[1][1][1])
								* (u[_i][_j + 1][_k] - u[_i][_j][_k])
								+ (_diff[1][0][1] + _diff[1][1][1])
								* (u[_i][_j - 1][_k] - u[_i][_j][_k])
								+ (_diff[1][1][2] + _diff[1][1][1])
								* (u[_i][_j][_k + 1] - u[_i][_j][_k]) + (_diff[1][1][0] + _diff[1][1][1])
								* (u[_i][_j][_k - 1] - u[_i][_j][_k]))
								* h2i + r;
						// update concentration (test for NaN)
						if (lop != lop) {
							System.out.println("---------------------");
							System.out
									.println("Multigrid problem (computing lop):");
							System.out.println("_g = " + _g);
							System.out.println("_i = " + _i);
							System.out.println("_j = " + _j);
							System.out.println("_k = " + _k);
							System.out.println("r = " + r);
							System.out.println("solute species:");
							Collection ss = Model.model().getSoluteSpecies();
							for (Iterator iter = ss.iterator(); iter.hasNext();) {
								SoluteSpecies s = (SoluteSpecies) iter.next();
								System.out.println(s._name + " = "
										+ s._mg[_g][_i][_j][_k]);
							}
							// particulates:
							System.out.println("particulate species:");
							Collection ps = Model.model()
									.getParticulateSpecies();
							for (Iterator iter = ps.iterator(); iter.hasNext();) {
								ParticulateSpecies p = (ParticulateSpecies) iter
										.next();
								System.out.println(p._name + " = "
										+ p._mg[_g][_i][_j][_k]);
							}
							System.out.println("---------------------");
							String str = "NaN generated in multigrid solver "
									+ "while computing rate for " + c._name;
							throw new ModelRuntimeException(str);
						}
						res._mg[_g][_i][_j][_k] = lop;
					}
			}
		}
		_boundaryConditions.refreshBoundaryConditions(res._mg[_g]);
	}

	/**
	 * Solves the reaction-diffusion equation by relaxation. Slower alternative
	 * to using method solve. NOTE: this method is to be used for testing
	 * purposes.
	 * 
	 * @param chem
	 * @param bac
	 */
	public static void solveByRelax(SoluteSpecies[] chem,
			ParticulateSpecies[] bac, BoundaryLayer bl)
			throws MultigridSystemNotSetException {
		_g = _order - 1;
		// create a relative diffusion data
		RelativeDiffusion relDiff = new RelativeDiffusion();
		bl.setBoundaryLayer(bac, _boundaryConditions);
		relDiff.computeValues(bac, bl, _boundaryConditions);
		// Initialize concentration of all solutes to value of
		// bulk concentration
		for (int i = 0; i < chem.length; i++) {
			chem[i].resetMultigridCopies();
		}
		// iterate for a guiven number of iterations
		for (int i = 0; i < 5000; i++)
			relax(chem, relDiff, bl);
	}

	/**
	 * Return values in finest grid as a matrix in a formatted string
	 * 
	 * @return string output
	 */
	public String finestGridToString() {
		return gridToString(_order - 1);
	}

	/**
	 * Return values in finest grid as a array of matrices in formatted strings
	 * 
	 * @return array of string output
	 */
	public String[] gridsToStrings() {
		String[] out = new String[_order];
		for (int i = 0; i < out.length; i++) {
			out[i] = gridToString(i);
		}
		return out;
	}

	/**
	 * Return values in grid g as a matrix in a formatted string
	 * 
	 * @return string output
	 */
	public String gridToString(int g) {
		return MultigridUtils.coreMatrixToString(_mg[g]);
	}

	/**
	 * @return the value of the finest grid
	 */
	public float[][][] getFinestGrid() {
		return _mg[_order - 1];
	}

	/**
	 * @return name of multigrid variable
	 */
	public String getName() {
		return _name;
	}

	/**
	 * Write a grid to file
	 * 
	 * @param g
	 *            grid index
	 * @param s
	 *            name of (path to) directory
	 * @throws IOException
	 *             if file may not be oppened
	 */
	protected void gridToFile(int g, String s) throws IOException {
		java.io.File f = new java.io.File(s + _name + g + ".txt");
		java.io.FileWriter fr = new java.io.FileWriter(f);
		fr.write(gridToString(g));
		fr.close();
	}

	/**
	 * Write out the finest grid to a file
	 * 
	 * @param s
	 *            directory name
	 * @throws IOException
	 */
	public void finestGridToFile(String s) throws IOException {
		gridToFile(_order - 1, s);
	}

	/**
	 * Write out every grid to files
	 * 
	 * @param s
	 *            directory name
	 * @throws IOException
	 */
	public void allGridsToFiles(String s) throws IOException {
		for (int i = 0; i < _order; i++)
			gridToFile(i, s);
	}

	/**
	 * Set the value of number of pre and post processing steps
	 * 
	 * @param npre
	 *            number of preprocessing steps
	 * @param npos
	 *            number of postprocessing steps
	 */
	public static void setSteps(int npre, int npos) {
		nPreSteps = npre;
		nPosSteps = npos;
	}

	/**
	 * Impose a profile for a solute species, such as in the case of the oxygen
	 * profile imposed in the Laspidou model
	 * 
	 * @param chem
	 * @param lambda
	 */
	public static void imposeProfile(SoluteSpecies chem, float lambda) {
		float[][][] profile = chem._mg[_order - 1];
		float cbulk = chem.getBulkConcentration();
		int n = profile.length;
		int m = profile[0].length;
		int l = profile[0][0].length;
		// find the boundary layer height
		int h = (int) (_model.getCurrentBiofilmHeight() / _voxelSide) + 1;
		for (int i = 0; i < n; i++) {
			// impose profile if
			float v = (i > h ? 1 : (float) (Math.exp(-lambda * (h - i + 0.5f)
					* _voxelSide)))
					* cbulk;
			for (int j = 0; j < m; j++)
				for (int k = 0; k < l; k++)
					profile[i][j][k] = v;
		}
	}

	/**
	 * Compute the biovolume from the multigrid representation of the
	 * particulate species
	 * 
	 * @param bac
	 *            array of particulate species
	 * @return the total volume occupied by the biofilm (biovolume)
	 */
	public static float computeBiovolume(ParticulateSpecies[] bac) {
		// count the voxels occupied by the biofilm
		int n = 0;
		// get the information concerning all the biomass
		for (int i = 1; i < _n + 1; i++)
			for (int j = 1; j < _m + 1; j++)
				for (int k = 1; k < _l + 1; k++)
					for (int sp = 0; sp < bac.length; sp++)
						if (bac[sp]._mg[_order - 1][i + 1][j + 1][k + 1] > 0) {
							n++;
							break;
						}
		// convert to volume
		return _voxelVolume * ((float) n);
	}

	/**
	 * @return Returns the volume of finest grid voxel
	 */
	public static float get_voxelVolume() {
		return _voxelVolume;
	}

	/**
	 * @param os
	 * @throws IOException
	 */
	public static void serializeStaticState(ObjectOutputStream os)
			throws IOException {
		os.writeInt(nPreSteps);
		os.writeInt(nPosSteps);
		os.writeInt(_l);
		os.writeInt(_m);
		os.writeInt(_n);
		os.writeInt(_order);
		os.writeFloat(_voxelVolume);
		os.writeFloat(_voxelSide);
		os.writeFloat(_referenceSystemSide);
		os.writeInt(_g);
		os.writeInt(_i);
		os.writeInt(_j);
		os.writeInt(_k);
		os.writeBoolean(_gridIsSet);
		os.writeObject(_boundaryConditions);
		os.writeObject(_model);
	}

	/**
	 * @param os
	 * @throws IOException
	 * @throws ClassNotFoundException
	 */
	public static void deserializeStaticState(ObjectInputStream os)
			throws IOException, ClassNotFoundException {
		nPreSteps = os.readInt();
		nPosSteps = os.readInt();
		_l = os.readInt();
		_m = os.readInt();
		_n = os.readInt();
		_order = os.readInt();
		_voxelVolume = os.readFloat();
		_voxelSide = os.readFloat();
		_referenceSystemSide = os.readFloat();
		_g = os.readInt();
		_i = os.readInt();
		_j = os.readInt();
		_k = os.readInt();
		_gridIsSet = os.readBoolean();
		_boundaryConditions = (BoundaryConditions) (os.readObject());
		_model = (Model) (os.readObject());
	}

	/**
	 * Sets the boundary conditions for this system
	 * 
	 * @param bc
	 *            The _boundaryConditions to set.
	 */
	public static void setBoundaryConditions(BoundaryConditions bc) {
		_boundaryConditions = bc;
	}

	/**
	 * @return a 3D matrix with the same size as finest grid
	 */
	public static boolean[][][] create3DBooleanMatrixWithFinnerResolution() {
		return new boolean[_n + 2][_m + 2][_l + 2];
	}

	/**
	 * an array with the size of the finer grid including padding
	 * 
	 * @return {_n+2, _m+2, _l+2}
	 */
	public static int[] getFinnerGridSize() {
		int[] c = { _n + 2, _m + 2, _l + 2 };
		return c;
	}
}
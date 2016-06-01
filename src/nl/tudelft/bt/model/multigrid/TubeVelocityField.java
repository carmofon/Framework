/*
 * Created on Jul 10, 2005
 */
package nl.tudelft.bt.model.multigrid;

import java.util.Collection;

import nl.tudelft.bt.model.Model;
import nl.tudelft.bt.model.exceptions.ModelRuntimeException;
import nl.tudelft.bt.model.exceptions.MultigridSystemNotSetException;
import nl.tudelft.bt.model.multigrid.boundary_conditions.BoundaryConditions;
import nl.tudelft.bt.model.multigrid.boundary_conditions.TubeBoundaryConditions;
import nl.tudelft.bt.model.multigrid.boundary_layers.BoundaryLayer;
import nl.tudelft.bt.model.profiles1d.Radial2DMultigridProfile;
import nl.tudelft.bt.model.util.ExtraMath;

/**
 * Implements the velocity field for tube reactors, where the velocity is
 * determined by solving a dimensionless Navier-Stokes equation. Dimensionless
 * velocity is then scaled so that the value for the volumetric flow rate is
 * obtained. The tube velocity field is a singleton.
 * 
 * @author jxavier
 */
public class TubeVelocityField extends MultigridVariable {
	private static TubeVelocityField _singleton;

	private float _flowRate;

	private float _viscosity;

	private float _tubeLength;

	private float _tubeRadius;

	private float _referenceDiffusivity;

	// a static 3D array for neighborhood computations in relax
	private static final float[][][] _diffVisc = new float[3][3][3];

	// right-hand-side for multigrid
	private MultigridVariable _rhs;

	// the shear stress (computed from gradient of velocity)
	private MultigridVariable _shearStress;

	// the shear stress (computed from gradient of velocity)
	private BoundaryLayer _boundaryLayer;

	private RelativeViscosity _relativeViscosity;

	private float _truncationError;

	private TubeBoundaryConditions.VelocityBoundaryConditions _velocityBoundaryConditions;

	private Radial2DMultigridProfile _velocityRadialProfile;
	private Radial2DMultigridProfile _shearRadialProfile;
	
	/**
	 * Creates a boundary layer which sets as interior of boundary layer any
	 * point where velocity is lower or equal to a defined threshold value
	 * (_velocityThresholdForBoundaryLayer)
	 * 
	 * @author xavier
	 */
	public class VelocityBoundaryLayer extends BoundaryLayer {
		/**
		 * @throws MultigridSystemNotSetException
		 */
		public VelocityBoundaryLayer() throws MultigridSystemNotSetException {
			super();
		}

		public void setBoundaryLayer(ParticulateSpecies[] b,
				BoundaryConditions bc) {
			//determine the threshold velocity
			float avereageVelocity = _flowRate
					/ ExtraMath.areaOfACircle(_tubeRadius);
			float lc = (float) Math.sqrt(_tubeLength * _referenceDiffusivity
					/ 4 / avereageVelocity);
			float thresholdVelocity = 2
					* avereageVelocity
					* (1 - ExtraMath.sq(_tubeRadius - lc)
							/ ExtraMath.sq(_tubeRadius));
			// compute the integral of the dimensionless velocity
			float[][][] u = _singleton._mg[_order - 1];
			float[][][] bl = _mg[_order - 1];
			float sum = 0;
			// set the values outside boundary layer to 1, 0 otherwise
			MultigridUtils.setValues(bl, 0);
			for (int i = 1; i < u.length - 1; i++)
				for (int j = 1; j < u[i].length - 1; j++)
					for (int k = 1; k < u[i][j].length - 1; k++)
						bl[i][j][k] = (u[i][j][k] > thresholdVelocity ? 1.0f
								: 0);
		}
	}

	/**
	 * @throws MultigridSystemNotSetException
	 */
	private TubeVelocityField() throws MultigridSystemNotSetException {
		//allocates multigrid variable
		super();
		//sets VelocityField attributes
		_name = "Velocity";
		_rhs = new MultigridVariable();
		_shearStress = new MultigridVariable();
		_shearStress._name = "ShearStress";
		_relativeViscosity = new RelativeViscosity();
		_boundaryLayer = new VelocityBoundaryLayer();
		_velocityRadialProfile = new Radial2DMultigridProfile(this);
		_shearRadialProfile = new Radial2DMultigridProfile(_shearStress);
	}

	/**
	 * @return the singleton instance of TubeVelocityField
	 */
	public static TubeVelocityField getInstance() {
		if (_singleton == null) {
			try {
				_singleton = new TubeVelocityField();
			} catch (MultigridSystemNotSetException e) {
				throw new ModelRuntimeException(
						"Trying to get instance of TubeVelocity "
								+ "before multigrid system is completely defined");
			}
		}
		return _singleton;
	}

	/**
	 * Solves the velocity field using Navier-Stokes equation and computed the
	 * shear stress from a gradient of the velocity
	 * 
	 * @throws MultigridSystemNotSetException
	 */
	public void solveVelocityFieldAndComputeShear()
			throws MultigridSystemNotSetException {
		solveDimensionlessNavierStokes();
		scaleTheVelocityValues();
		computeShearStress();
		// update the values of the profiles
		_velocityRadialProfile.computeProfile();
		_shearRadialProfile.computeProfile();
	}

	/**
	 * Scale the values of the dimensionless velocity so that the integral of
	 * the velocity iver the area where flow is passing equals the volumetric
	 * flow rate. This operation is carried out ounly for the finnest grid.
	 * Values at courser grid levels are not scale (remain dimensionless) as
	 * they are not used for further computations
	 */
	private void scaleTheVelocityValues() {
		// compute the integral of the dimensionless velocity
		float[][][] u = _mg[_order - 1];
		float sum = 0;
		for (int i = 1; i < u.length - 1; i++)
			for (int j = 1; j < u[i].length - 1; j++)
				for (int k = 1; k < u[i][j].length - 1; k++)
					sum += u[i][j][k] * ExtraMath.sq(_voxelSide);
		//compute the scalling factor
		float f = _flowRate / sum;
		// scale the velocity
		for (int i = 1; i < u.length - 1; i++)
			for (int j = 1; j < u[i].length - 1; j++)
				for (int k = 1; k < u[i][j].length - 1; k++)
					u[i][j][k] *= f;
	}

	/**
	 * Computes the local shear stress at finner grid from the gradient of the
	 * velocity
	 */
	private void computeShearStress() {
		//get the grid node size
		float stepSize = _voxelSide;
		// compute the integral of the dimensionless velocity
		float[][][] u = _mg[_order - 1];
		float[][][] shear = _shearStress._mg[_order - 1];
		float sum = 0;
		for (int i = 1; i < u.length - 1; i++)
			for (int j = 1; j < u[i].length - 1; j++)
				for (int k = 1; k < u[i][j].length - 1; k++) {
					//TODO testing phase local versus water velocity
					//float localViscosity = _viscosity
					//* _relativeViscosity._mg[_g][_i][_j][_k];
					float localViscosity = _viscosity;
					// x-direction
					float dx = ExtraMath.maxSquare(u[i + 1][j][k] - u[i][j][k],
							u[i][j][k] - u[i - 1][j][k]);
					// y-direction
					float dy = ExtraMath.maxSquare(u[i][j + 1][k] - u[i][j][k],
							u[i][j][k] - u[i][j - 1][k]);
					// the norm of the shear vector
					shear[i][j][k] = (float) Math.sqrt(dx + dy) / _voxelSide
							* localViscosity;
				}
	}

	/**
	 * Solve the velocity flow in tube by solving Navier-Stokes equation in
	 * tube.
	 * 
	 * @throws MultigridSystemNotSetException
	 *             when system is not set
	 */
	private void solveDimensionlessNavierStokes()
			throws MultigridSystemNotSetException {
		//instanciate the boundary conditions
		_velocityBoundaryConditions = ((TubeBoundaryConditions) _boundaryConditions)
				.getVelocityBoundaryConditions();
		//get the references to the biomass species so that boundary conditions
		//may be obtained
		Model.model().updateBioDiscreteData();
		Collection particulates = Model.model().getParticulateSpecies();
		// two temporary multigrid variables are needed for the
		// computation, are initialized with 0 allover
		MultigridVariable itemp = new MultigridVariable();
		MultigridVariable itau = new MultigridVariable();
		// create a relative viscosity matrix:
		// 1 in liquid
		// very high in biofilm and carrier (to simulate 0 velocity flux)
		_relativeViscosity.computeValues(particulates,
				_velocityBoundaryConditions);
		_relativeViscosity.updateMultigridCopies();
		// Initialize velocity field
		resetMultigridCopies();
		// solve chemical concentrations on coarsest grid
		solveCoarsest();
		// nested iteration loop
		for (int outer = 1; outer < _order; outer++) {
			_g = outer;
			MultigridUtils.interpolate(_mg[_g], _mg[_g - 1],
					_velocityBoundaryConditions);
			// set each chemical's r.h.s. to 0
			MultigridUtils.setValues(_rhs._mg[_g], 0.0f);
			float iterror = 1.0f;
			// V-cycle loop
			for (int v = 0; v < VCYCLES; v++) {
				// downward stroke of V
				while (_g > 0) {
					// pre-smoothing
					for (int j = 0; j < nPreSteps; j++) {
						relax();
					}
					// restrict uh to uH
					MultigridUtils.restrict(_mg[_g], _mg[_g - 1],
							_velocityBoundaryConditions);
					//
					lop(itemp);
					//
					MultigridUtils.restrict(itemp._mg[_g], itemp._mg[_g - 1],
							_velocityBoundaryConditions);
					// reduce grid value _g temporarily
					_g--;
					lop(itau);
					MultigridUtils.subtractTo(itau._mg[_g], itemp._mg[_g]);
					// sum tau to rhs of _g - 1
					MultigridUtils.restrict(_rhs._mg[_g + 1], _rhs._mg[_g],
							_velocityBoundaryConditions);
					MultigridUtils.addTo(_rhs._mg[_g], itau._mg[_g]);
					// compute the truncation error for this V-cycle
					// for all chemicals
					if (_g + 1 == outer)
						_truncationError = ALPHA
								* MultigridUtils.computeNorm(itau._mg[_g]);
				}
				// bottom of V
				solveCoarsest();
				// upward stroke of V
				while (_g < outer) {
					_g++;
					MultigridUtils.restrict(_mg[_g], itemp._mg[_g - 1],
							_velocityBoundaryConditions);
					MultigridUtils.subtractTo(_mg[_g - 1], itemp._mg[_g - 1]);
					MultigridUtils.interpolate(itau._mg[_g], _mg[_g - 1],
							_velocityBoundaryConditions);
					MultigridUtils.addTo(_mg[_g], itau._mg[_g]);
					// post-smoothing
					for (int j = 0; j < nPosSteps; j++) {
						relax();
					}
				}
				// break the V-cycles if remaining error is dominated
				// by local truncation error (see p. 884 of Numerical Recipes)
				boolean breakVCycle = true;
				// compute the residue for this solute species
				lop(itemp);
				MultigridUtils.subtractTo(itemp._mg[_g], _rhs._mg[_g]);
				float res = MultigridUtils.computeNorm(itemp._mg[_g]);
				// confirm that criterium is met
				if (v > 0) {
					System.out.println("grid " + _g + "; v " + v + "; " + _name
							+ " res " + res + "; truncerr " + _truncationError);
				}
				// confirm that criterium is met for each solute
				if (res > _truncationError) {
					breakVCycle = false;
					break;
				}
				if (breakVCycle)
					break;
			}
		}

	}

	/**
	 * Find solution for the coarsest grid. Sets the current grid to coarsest,
	 * velocity values to 0 and relaxes NSOLVE times.
	 */
	private void solveCoarsest() {
		_g = COARSEST;
		float err = 0;
		// reset coarsest grid to bulk concentration
		MultigridUtils.setValues(_mg[COARSEST], 0);
		// relax NSOLVE times
		for (int j = 0; j < NSOLVE; j++) {
			relax();
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
	private void relax() {
		// stores the number of elments that
		// do not comply to precision criterium
		int elementsToGo = 0;
		float r, dr;
		int n = _relativeViscosity._mg[_g].length - 2;
		int m = _relativeViscosity._mg[_g][0].length - 2;
		int l = _relativeViscosity._mg[_g][0][0].length - 2;
		// set sizes of the system (dimensionless)
		float h = 1 / ((float) n - 1);
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
						// Case: Inside boundary layer
						// Equations must be solved here
						float[][][] u = _mg[_g];
						r = 1.0f;
						dr = 0.0f;
						// compute diffusivity values
						// and that of surrounding neighbors
						float dc = 1.0f;
						_diffVisc[0][1][1] = _relativeViscosity._mg[_g][_i - 1][_j][_k];
						_diffVisc[2][1][1] = _relativeViscosity._mg[_g][_i + 1][_j][_k];
						_diffVisc[1][0][1] = _relativeViscosity._mg[_g][_i][_j - 1][_k];
						_diffVisc[1][2][1] = _relativeViscosity._mg[_g][_i][_j + 1][_k];
						_diffVisc[1][1][0] = _relativeViscosity._mg[_g][_i][_j][_k - 1];
						_diffVisc[1][1][2] = _relativeViscosity._mg[_g][_i][_j][_k + 1];
						_diffVisc[1][1][1] = _relativeViscosity._mg[_g][_i][_j][_k];
						// compute L operator
						float lop = ((_diffVisc[2][1][1] + _diffVisc[1][1][1])
								* (u[_i + 1][_j][_k] - u[_i][_j][_k])
								+ (_diffVisc[0][1][1] + _diffVisc[1][1][1])
								* (u[_i - 1][_j][_k] - u[_i][_j][_k])
								+ (_diffVisc[1][2][1] + _diffVisc[1][1][1])
								* (u[_i][_j + 1][_k] - u[_i][_j][_k])
								+ (_diffVisc[1][0][1] + _diffVisc[1][1][1])
								* (u[_i][_j - 1][_k] - u[_i][_j][_k])
								+ (_diffVisc[1][1][2] + _diffVisc[1][1][1])
								* (u[_i][_j][_k + 1] - u[_i][_j][_k]) + (_diffVisc[1][1][0] + _diffVisc[1][1][1])
								* (u[_i][_j][_k - 1] - u[_i][_j][_k]))
								* h2i + r;
						// compute derivative of L operator
						float dlop = -h2i
								* (6.0f * _diffVisc[1][1][1]
										+ _diffVisc[2][1][1]
										+ _diffVisc[0][1][1]
										+ _diffVisc[1][2][1]
										+ _diffVisc[1][0][1]
										+ _diffVisc[1][1][2] + _diffVisc[1][1][0])
								+ dr;
						// compute residual
						float res = (lop - _rhs._mg[_g][_i][_j][_k]) / dlop;
						// update concentration (test for NaN)
						if (res != res) {
							throw new ModelRuntimeException(
									"NaN generated in multigrid solver"
											+ "while computing rate for"
											+ _name);
						}
						u[_i][_j][_k] -= res;
						//if negative concentrations, put 0 value
						u[_i][_j][_k] = (u[_i][_j][_k] < 0 ? 0 : u[_i][_j][_k]);
					}
				}
			}
			// refresh the padding elements to enforce
			// boundary conditions for all solutes
			_velocityBoundaryConditions.refreshBoundaryConditions(_mg[_g]);
		}
	}

	/**
	 * Compute the L-operator
	 * 
	 * @param res
	 *            residual
	 * @param visc
	 *            relative viscosity
	 */
	private void lop(MultigridVariable res) {
		int n = _relativeViscosity._mg[_g].length - 2;
		int m = _relativeViscosity._mg[_g][0].length - 2;
		int l_ = _relativeViscosity._mg[_g][0][0].length - 2;
		float h = 1.0f / ((float) n - 1);
		float h2i = 0.5f / (h * h);
		// iterate through system
		for (_k = 1; _k <= l_; _k++) {
			for (_j = 1; _j <= m; _j++) {
				for (_i = 1; _i <= n; _i++) {
					// for simplification and easier access to
					// the current grid:
					float[][][] u = _mg[_g];
					// current rate for this solute
					float r = 1.0f;
					// compute diffusivity values
					// and that of surrounding neighbors
					float dc = 1.0f;
					_diffVisc[0][1][1] = _relativeViscosity._mg[_g][_i - 1][_j][_k];
					_diffVisc[2][1][1] = _relativeViscosity._mg[_g][_i + 1][_j][_k];
					_diffVisc[1][0][1] = _relativeViscosity._mg[_g][_i][_j - 1][_k];
					_diffVisc[1][2][1] = _relativeViscosity._mg[_g][_i][_j + 1][_k];
					_diffVisc[1][1][0] = _relativeViscosity._mg[_g][_i][_j][_k - 1];
					_diffVisc[1][1][2] = _relativeViscosity._mg[_g][_i][_j][_k + 1];
					_diffVisc[1][1][1] = _relativeViscosity._mg[_g][_i][_j][_k];
					// compute L operator
					res._mg[_g][_i][_j][_k] = ((_diffVisc[2][1][1] + _diffVisc[1][1][1])
							* (u[_i + 1][_j][_k] - u[_i][_j][_k])
							+ (_diffVisc[0][1][1] + _diffVisc[1][1][1])
							* (u[_i - 1][_j][_k] - u[_i][_j][_k])
							+ (_diffVisc[1][2][1] + _diffVisc[1][1][1])
							* (u[_i][_j + 1][_k] - u[_i][_j][_k])
							+ (_diffVisc[1][0][1] + _diffVisc[1][1][1])
							* (u[_i][_j - 1][_k] - u[_i][_j][_k])
							+ (_diffVisc[1][1][2] + _diffVisc[1][1][1])
							* (u[_i][_j][_k + 1] - u[_i][_j][_k]) + (_diffVisc[1][1][0] + _diffVisc[1][1][1])
							* (u[_i][_j][_k - 1] - u[_i][_j][_k]))
							* h2i + r;
				}
			}
		}
		_velocityBoundaryConditions.refreshBoundaryConditions(res._mg[_g]);
	}

	/**
	 * Initialize the values of the velocity to 0 at all grid levels
	 */
	private void resetMultigridCopies() {
		for (int i = _order - 1; i > COARSEST; i--) {
			MultigridUtils.setValues(_mg[i], 0);
		}
	}

	/**
	 * Set the value of the volumetric flow rate
	 * 
	 * @param rate
	 *            The _flowRate to set.
	 * @param viscosity
	 *            the liquid's viscosity
	 * @param tubeLength
	 *            TODO
	 */
	public void setFlowProperties(float rate, float viscosity, float tubeRadius,
			float tubeLength, float referenceDiffusivity) {
		_flowRate = rate;
		_viscosity = viscosity;
		_tubeLength = tubeLength;
		_tubeRadius = tubeRadius;
		_referenceDiffusivity = referenceDiffusivity;
	}

	/**
	 * @return Returns the _shearStress.
	 */
	public MultigridVariable getShearStress() {
		return _shearStress;
	}

	/**
	 * @return Returns the _boundaryLayer.
	 */
	public BoundaryLayer getBoundaryLayer() {
		return _boundaryLayer;
	}
	/**
	 * @return Returns the instance of the velocity radial profile.
	 */
	public Radial2DMultigridProfile getVelocityRadialProfile() {
		return _velocityRadialProfile;
	}
	/**
	 * @return Returns the _shearRadialProfile.
	 */
	public Radial2DMultigridProfile getShearRadialProfile() {
		return _shearRadialProfile;
	}
}
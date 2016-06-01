package nl.tudelft.bt.model.work.will;

import nl.tudelft.bt.model.multigrid.SoluteSpecies;
import nl.tudelft.bt.model.ContinuousCoordinate;
import nl.tudelft.bt.model.Model;
import nl.tudelft.bt.model.particlebased.BiomassParticle;
import nl.tudelft.bt.model.reaction.ProcessFactor;
import nl.tudelft.bt.model.reaction.Saturation;
import nl.tudelft.bt.model.util.ExtraMath;
import nl.tudelft.bt.model.util.ContinuousCoordinateOps;

/** Implements chemotaxis/directed migration in a biomass particle
 * 
 * @author Will Chang - 2010
 *
 */

public class ChemotacticBiomassParticle extends BiomassParticle {
	
	protected ChemotacticBiomassSpecies _motileSpecies;
	
	protected SoluteSpecies _attractantSpecies;
	
	protected float _attractantAtCenter;
	
	// the velocity of the particle in the previous iteration
	protected ContinuousCoordinate _vp = new ContinuousCoordinate();
	
//	protected float _maximumTimeStep;
	
	// checkers for whether cell is at a positive grid boundary
	protected boolean _atXEdge;
	
	protected boolean _atYEdge;
	
	protected boolean _atZEdge;
	
	protected ProcessFactor _motilityFactor;
	
	/**
	 * Constructor
	 * 
	 * @param s the chemotactic biomass species
	 */
	public ChemotacticBiomassParticle(ChemotacticBiomassSpecies s) {
		super(s);
		_motileSpecies = s;
		_attractantSpecies = s.getChemoattractant();
		_attractantAtCenter = _attractantSpecies.getValueAt(getCenter());
		_motilityFactor = s._motilityFactor;
	}
	
	@Override
	public float grow(float t) {
		// get the cell center coordinates before growth and movement
		ContinuousCoordinate cp = new ContinuousCoordinate(_centerX, _centerY, _centerZ);
		cp.scale(-1);
		// grow and move, shove etc.
		ContinuousCoordinate movementStep = getMovementStep(t);
		move(movementStep.x, movementStep.y, movementStep.z);
		float m = super.grow(t);
		// get the new cell center and retroactively calculate velocity
		// update the velocity
		// don't perform on first iteration to avoid artifacts from spreading
		// after inoculating
		if (Model.model().getIterationCounter() > 1) {
			ContinuousCoordinate cn = getCenter();
			ContinuousCoordinate d = ContinuousCoordinateOps.sum(
					new ContinuousCoordinate[]{ cn, cp }
			);
			// update the velocity and adhesive radius
			_vp = ContinuousCoordinateOps.newScale(d, 1/t);
		}
		return m;
	}
	
	@Override
	/**
	 * Code dupe from GranuleBiomassParticle - maybe clean up later
	 * Set position of center. All borders are cyclic
	 * 
	 * @param x
	 *            x [micron]
	 * @param y
	 *            y [micron]
	 * @param z
	 *            z [micron]
	 */
	public void setCenter(float x, float y, float z) {
		// Implement periodic boundaries
		float Lxx = _model.systemSize.x;
		if (x >= Lxx)
			_centerX = x - Lxx;
		else if (x < 0)
			_centerX = Lxx + x;
		else
			_centerX = x;
		float Lyy = _model.systemSize.y;
		if (y >= Lyy)
			_centerY = y - Lyy;
		else if (y < 0)
			_centerY = Lyy + y;
		else
			_centerY = y;
		float Lzz = _model.systemSize.z;
		if (Lzz == 0)
			// 2D case
			_centerZ = 0;
		else if (z >= Lzz)
			_centerZ = z - Lzz;
		else if (z < 0)
			_centerZ = Lzz + z;
		else
			_centerZ = z;
		setShovingMapPositionAndAddToShovingMap();
	}
	
	/**
	 * Calculates chemotaxis movement vector in one iteration of time length t
	 * 
	 * @param t time length
	 */
	protected ContinuousCoordinate getChemotaxisDirection() {
		// get the coordinates of cell center
		_attractantAtCenter = _attractantSpecies.getValueAt(getCenter());
		// get the gradient of the chemoattractant
		ContinuousCoordinate gradient = computeChemoattractantGradient();
		float gradMagnitude = ContinuousCoordinateOps.getLength(gradient);
		ContinuousCoordinate dir;
		// choose between the correct and random directions
		if (chooseCorrectDirection(gradMagnitude)) {
			gradient.normalize();
			dir = gradient;
		} else {
			// generate the alternative random direction vector
			dir = ContinuousCoordinateOps.getRandomUnitVector();
		}
		return dir;
	}
	
	/**
	 * Computes the scalar length of movement in iteration of length t
	 * 
	 * @param t
	 * @return
	 */
	protected float getChemotaxisLength(float t) {
		float d = _motileSpecies.getDiffusivity();
		// without Gaussian scaling
//		float factor = _motilityFactor.getValue();
		// with Gaussian scaling
		float factor = _motilityFactor.getValue() * 
			(Model.model().getRandomFromNormalDistribution() + 1f); 
		// with power law scaling
//		float factor = computePowerLawSpeedFactor();
		float l = ExtraMath.sqrt(d * t) * factor;
		return -l;
	}
	
	/**
	 * Computes the absolute movement vector
	 * 
	 * @param t iteration length
	 * @return vector of movement
	 */
	protected ContinuousCoordinate getMovementStep(float t) {
		// code for motility
		ContinuousCoordinate dir = getChemotaxisDirection();
		float length = getChemotaxisLength(t);
		ContinuousCoordinate movementStep = ContinuousCoordinateOps.newScale(dir, length);
		return movementStep;
	}
		
	/**
	 * First-order gradient calculation
	 * @return
	 * 		float[] the absolute gradient of the attractant species 
	 */
	protected ContinuousCoordinate computeChemoattractantGradient() {
		float voxelSide = Model.model().getGridElementSide();
		// get the attractant concentration in this voxel
		float thisConcentration = _attractantSpecies.getValueAt(getCenter());
		float delX;
		float delY;
		float delZ;
		// check whether cell is at a positive edge of grid
		checkEdges();
		if (!_atXEdge) {
			delX = voxelSide;
		} else {
			delX = -voxelSide;
		}
		ContinuousCoordinate cX = new ContinuousCoordinate ((_centerX + delX), _centerY, _centerZ);
		float gradX = (_attractantSpecies.getValueAt(cX) - thisConcentration) / delX;
		if (!_atYEdge) {
			delY = voxelSide;
		} else {
			delY = -voxelSide;
		}
		ContinuousCoordinate cY = new ContinuousCoordinate (_centerX, (_centerY + delY), _centerZ);
		float gradY = (_attractantSpecies.getValueAt(cY) - thisConcentration) / delY;
		float gradZ;
		if (Model.model().getDimensionality() == 3) {
			if (!_atZEdge) {
				delZ = voxelSide;
			} else {
				delZ = -voxelSide;
			}
			ContinuousCoordinate cZ = new ContinuousCoordinate (_centerX, _centerY, (_centerZ + delZ));
			gradZ = (_attractantSpecies.getValueAt(cZ) - thisConcentration) / delZ;
			gradZ = gradZ / delZ;
		} else {gradZ = 0;};
		ContinuousCoordinate gradient = new ContinuousCoordinate(gradX, gradY, gradZ);
		return gradient;
	}
	
	/**
	 * Decides whether to take the "right" direction for migration.
	 * 
	 * @param gradMag
	 * 		the magnitude of the gradient vector
	 * @param conc
	 * 		the concentration of chemoattractant at cell center
	 * @return
	 * 		TRUE if take the right direction this iteration, FALSE if
	 * 		take a random direction
	 */
	protected boolean chooseCorrectDirection(float gradMag) {
		if (_attractantAtCenter < _motileSpecies.getKspeed()) {
			return false;
		}
		// fractional gradient across cell diameter = 
		// gradient / concentration * cell diameter
		float fG = (gradMag / _attractantAtCenter) * 2 * this.getRadius();
		// calculate the probability of choosing the correct orientation
//		float p = fG / (fG + _motileSpecies.getKspread());
		double p = Math.pow(fG, 2) / (Math.pow(fG, 2) + Math.pow(_motileSpecies.getKspread(), 2));
		// compare the probability to a uniform random number
		// choose the correct orientation if the probability is greater
		boolean chooseCorrect = ((float) p) > Model.model().getRandom();
		return chooseCorrect;
	}
	
	/**
	 * Updates whether the cell is at the positive X, Y, or Z edges
	 */
	protected void checkEdges() {
		float voxelSide = Model.model().getGridElementSide();
		if (_centerX + voxelSide > Model.model().systemSize.x) {
			_atXEdge = true;
		} else {_atXEdge = false;}
		if (_centerY + voxelSide > Model.model().systemSize.y) {
			_atYEdge = true;
		} else {_atYEdge = false;}
		if (_centerZ + voxelSide > Model.model().systemSize.z) {
			_atZEdge = true;
		} else {_atZEdge = false;}
	}
	
	/**
	 * @return the velocity of this cell in the previous turn
	 */
	public ContinuousCoordinate getVelocity() {
		return _vp;
	}
	
	/**
	 * Override the time step calculation to account for migration
	 * (do not migrate more than 1 cell radius in a time step)
	 * 
	 * @return the growth time step OR migration time step, whichever one is
	 * smaller
	 */
	@Override
	public float getMaximumTimeStep() {
		float growthTimeStep = super.getMaximumTimeStep();
		float dispTime = getMigrationTimeStep();
		return Math.min(growthTimeStep, dispTime);
	}
	
	/**
	 * Calculate the maximum time step without displacing the particle more
	 * than its own radius (assuming constant velocity)
	 * 
	 * @return time to displace particle by r
	 */
	public float getMigrationTimeStep() {
		float dispTime;
		// compute the individual displacement time on the first iteration
		if (Model.model().getIterationCounter() <= 1) {
			float factor = _motilityFactor.getValue();
			double dd = Math.pow(_r, 2) / 
				(_motileSpecies.getDiffusivity() * Math.pow(factor, 2));
			dispTime = (float) dd;
		}
		else {
			dispTime = _r / ContinuousCoordinateOps.getLength(_vp);
		}
		/* 
		 if velocity is 0 dispTime will be +infinity, there will be an error
		 if there is also no biomass growth
		 handle this exception by forcing the max time step of the model
		 (this must be set to > 0!)
		 * 
		 */
		if (dispTime == Float.POSITIVE_INFINITY) {
			dispTime = Model.model().getMaximumTimeStep();
		}
		return dispTime;
	}
}

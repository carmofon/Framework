package nl.tudelft.bt.model.work.benzoate;

import java.awt.Color;

import nl.tudelft.bt.model.BiomassSpecies;
import nl.tudelft.bt.model.Model;
import nl.tudelft.bt.model.particlebased.BiomassParticle;
import nl.tudelft.bt.model.util.ExtraMath;

/**
 * Implements surface motility into a biomass particle
 * 
 * @author Joao Xavier (jxavier@cgr.harvard.edu) - Oct 12, 2006
 */
public class MotileBiomassParticle extends BiomassParticle {
	private MotileBiomassSpecies _motileSpecies;

	private boolean _isMotile;

	public MotileBiomassParticle(MotileBiomassSpecies s) {
		super(s);
		_motileSpecies = s;
	}

	@Override
	public float grow(float t) {
		float m = super.grow(t);
		// code for motility
		moveBySurfaceMotility(t);
		return m;
	}

	/**
	 * Implements the motility
	 * 
	 * @param t
	 *            the time step of present iteration
	 */
	private void moveBySurfaceMotility(float t) {
		// determine a random move direction (2D, since movement is on surface)
		float moveTheta = Model.model().getRandom() * ExtraMath.PI2;
		float d = _motileSpecies.getDiffusivity();
		// following lines reduce themotility diffusivity as a funciton of the
		// growth
		// float d = _motileSpecies.getDiffusivity()
		// * (1 - _composition.getRelativeGrowth());
		// determine movement length from normally distributed random variable
		float length = Model.model().getRandomFromNormalDistribution()
				* ExtraMath.sqrt(d * t);
		// create movement vector
		float delx = 0; // is the vertical axis
		float dely = length * (float) (Math.cos(moveTheta));
		float delz = length * (float) (Math.sin(moveTheta));
		// check if bacteria is stuck
		// by checking if there are neighbors above
		if (!hasNeighborOnTop()) {
			// move
			move(0, -dely, -delz);
			delx = distanceToNearestNeighborBelow();
			move(delx, 0, 0);
			delx = distanceToNearestEmptyPlaceAbove();
			move(-delx, 0, 0);
			_isMotile = true;
		} else {
			delx = distanceToNearestNeighborBelow();
			move(delx, 0, 0);
			_isMotile = false;
		}
	}

	@Override
	public Color getColorCore() {
		if (_isMotile)
			return Color.blue;
		return super.getColorCore();
	}

	@Override
	protected Object clone() throws CloneNotSupportedException {
		// TODO Auto-generated method stub
		return super.clone();
	}
}

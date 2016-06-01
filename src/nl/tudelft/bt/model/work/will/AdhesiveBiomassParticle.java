package nl.tudelft.bt.model.work.will;

import java.util.ArrayList;

import nl.tudelft.bt.model.exceptions.NonMatchingNumberException;
import nl.tudelft.bt.model.particlebased.BiomassParticle;
import nl.tudelft.bt.model.particlebased.BiomassParticleSet;
import nl.tudelft.bt.model.ContinuousCoordinate;
import nl.tudelft.bt.model.Model;
import nl.tudelft.bt.model.util.ContinuousCoordinateOps;
import nl.tudelft.bt.model.util.ExtraMath;

public class AdhesiveBiomassParticle extends ChemotacticBiomassParticle {
	
	protected AdhesiveBiomassSpecies _adhesiveSpecies;
		
	private float _ws;
	
	private float _wg;
	
//	// the cell movement velocity vector of the previous turn
//	private ContinuousCoordinate _vp = new ContinuousCoordinate();
//	
//	// the cell movement velocity vector of the current turn
//	private ContinuousCoordinate _vn = new ContinuousCoordinate();
	
	// the relative range of cell adhesion
	private float _ra;
	
	// the absolute range of cell adhesion
	private float _adhesiveRadius;
	
	// the adhesion target cells
	private BiomassParticleSet _adhesionTargets = new BiomassParticleSet();
		
	/**
	 * Implements cell-cell adhesion in a chemotactic cell
	 * 
	 * @param s
	 * 	the adhesive biomass species
	 * 
	 * TODO: somehow make it possible to specify which cell types to adhere to
	 * */
	
	public AdhesiveBiomassParticle(AdhesiveBiomassSpecies s) {
		super(s);
		_adhesiveSpecies = s;
		_ra = s.getra();
		setAdhesiveRadius();
		_ws = _adhesiveSpecies.getws();
		_wg = _adhesiveSpecies.getwg();
	}
	
	/**
	 * Resets the list of adhesion targets, then
	 * goes through the cells in neighboring voxels of the shoving grid,
	 * if they are AdhesiveBiomassParticles, check to see if they are
	 * within adhesion range (= pairwise sum of adhesive radii),
	 * if they are, add them to the list of adhesion targets. 
	 */
	public void getAdhesionTargets() {
		_adhesionTargets.reset();
		BiomassParticleSet adhesionNeibs = 
			_particleContainer.getPotentialShovers(this);
		while (adhesionNeibs.hasNext()) {
			BiomassParticle neib = adhesionNeibs.next();
			// currently all adhesive cells stick to all others
			if ((neib instanceof AdhesiveBiomassParticle) & (neib != this)) {
				float range = 
					((AdhesiveBiomassParticle) neib).getAdhesiveRadius()
					+ this._adhesiveRadius;
				float d = ExtraMath.pointDistance(neib.getCenter(), 
						this.getCenter());
				// if the distance between two cell centers is smaller than 
				// the sum of adhesive radii, add the neighbor cell to the
				// list of adhesion targets
				if (d < range) {
					_adhesionTargets.add(neib);
				}
			}
		}
	}
	
	/**
	 * Calculates the contribution of a neighboring cell's position and velocity
	 * to focal cell's movement this turn
	 * 
	 * @param neib a neighboring AdhesiveBiomassParticle on the shoving grid
	 * @param t time step of current iteration
	 * @return the net contribution vector from the neighboring cell
	 */
	public ContinuousCoordinate getAdhesionContribution(AdhesiveBiomassParticle neib, float t) {
		// the velocity of the neighbor
		ContinuousCoordinate nc = neib.getCenter();
		ContinuousCoordinate cji = new ContinuousCoordinate(nc.x - _centerX, nc.y - _centerY, nc.z = _centerZ);
		float s = ContinuousCoordinateOps.getLength(cji) - (neib.getRadius() + this._r);
		cji.normalize();
		cji.scale(-s);
		ContinuousCoordinate vj = neib.getVelocity();
		ContinuousCoordinate rj = new ContinuousCoordinate(vj.x - _vp.x, vj.y - _vp.y, vj.z - _vp.z);
		ContinuousCoordinate dj = ContinuousCoordinateOps.sum(new ContinuousCoordinate[] { 
				ContinuousCoordinateOps.newScale(rj, t), cji });
		return dj;
	}
	
	/**
	 * Return the average adhesion contribution from neighbors
	 * (Should only be called when any neighbors present)
	 * 
	 * @param neibs set of neighboring AdhesiveBiomassParticles
	 * @param t time step of current iteration
	 * @return normalized adhesion direction
	 * @throw NonMatchingNumberException 
	 * 			if number of neighbors <= 0
	 */
	public ContinuousCoordinate getTotalAdhesion(BiomassParticleSet neibs, float t) 
		throws NonMatchingNumberException {
		int nneibs = neibs.getNumberOfParticles();
		// prevent divide by 0
		if (nneibs == 0) {
			throw new NonMatchingNumberException(
					"Tried to get adhesion of 0 neighbors!");
		}
		ContinuousCoordinate ds = new ContinuousCoordinate();
		while (neibs.hasNext()) {
			AdhesiveBiomassParticle neib = (AdhesiveBiomassParticle) neibs.next();
			ContinuousCoordinate dj = getAdhesionContribution(neib, t);
			ds.add(dj);
		}
		ds.scale(1/nneibs);
		return ds;
	}
	
	/**
	 * Computes the individual and social contributions to migration direction
	 * and returns the weighted sum of the contributions as a
	 * unit-length ContinuousCoordinate
	 * 
	 * Overrides the chemotaxis-only method in ChemotacticBiomassParticle
	 * with a combined chemotaxis-adhesion method
	 * 
	 * @param t time step of current iteration
	 * @return ContinuousCoordinate of unit length
	 */
	@Override
	protected ContinuousCoordinate getMovementStep(float t) {
		ContinuousCoordinate dg = ContinuousCoordinateOps.newScale(getChemotaxisDirection(),
				getChemotaxisLength(t));
		getAdhesionTargets();
		int nneibs = _adhesionTargets.getNumberOfParticles();
		// calculate the average contribution to migration from velocities of neighbors
		// (if any)
		ContinuousCoordinate ds = new ContinuousCoordinate();
		if (nneibs > 0) {
			ds = getTotalAdhesion(_adhesionTargets, t);
		}
		// weighted sum of individual and group contributions to migration
		ds.scale(_ws);
		dg.scale(_wg);
		ContinuousCoordinate movementStep = ContinuousCoordinateOps.sum(
				new ContinuousCoordinate[] { ds, dg }
		);
		return movementStep;
	}
	
	@Override
	public float grow(float t) {
//		// get the cell center coordinates before growth and movement
//		ContinuousCoordinate cp = new ContinuousCoordinate(_centerX, _centerY, _centerZ);
//		cp.scale(-1);
//		// grow and move, shove etc.
//		// movement methods are called in super method
		float m = super.grow(t);
//		// get the new cell center and retroactively calculate velocity
//		// update the velocity
//		// don't perform on first iteration to avoid artifacts from spreading
//		// after inoculating
//		if (Model.model().getIterationCounter() > 1) {
//			ContinuousCoordinate cn = getCenter();
//			ContinuousCoordinate d = ContinuousCoordinateOps.sum(
//					new ContinuousCoordinate[]{ cn, cp }
//			);
//			// update the velocity and adhesive radius
//			_vp = ContinuousCoordinateOps.newScale(d, 1/t);
//		}
		setAdhesiveRadius();
		return m;
	}
	
	/**
	 * Resets the adhesive radius given the current cell radius
	 * (do this after growth)
	 */
	protected void setAdhesiveRadius() {
		_adhesiveRadius = _ra * this.getRadius();
	}
	
	/**
	 * Gets the effective adhesive radius of this cell
	 * 
	 * @return _adhesiveRadius
	 */
	public float getAdhesiveRadius() {
		return _adhesiveRadius;
	}
	
//	/**
//	 * @return the velocity of this cell in the previous turn
//	 */
//	public ContinuousCoordinate getVelocity() {
//		return _vp;
//	}

}

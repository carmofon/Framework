/**
 * 
 */
package nl.tudelft.bt.model.work.will;

import java.awt.Color;

import nl.tudelft.bt.model.ContinuousCoordinate;
import nl.tudelft.bt.model.Model;
import nl.tudelft.bt.model.particlebased.BiomassParticle;
import nl.tudelft.bt.model.particlebased.BiomassParticleSet;
import nl.tudelft.bt.model.util.ContinuousCoordinateOps;
import nl.tudelft.bt.model.util.ExtraMath;

/**
 * @author wkc
 *
 */
public class ZonalVicsekBiomassParticle extends MutatorBiomassParticle {
	
	private ZonalVicsekBiomassSpecies _species;

	// stop simulation on periodic boundary crossing?
	private boolean _stopOnPeriodicBoundaryCrossing;
	
	// migrating this turn?
	protected boolean _migratingThisIteration = false;
	
	// the polarization vector
	private ContinuousCoordinate _polarization;
	
	// gradient sensing
	private ContinuousCoordinate _pg = new ContinuousCoordinate();
	
	// neighbors
	private ContinuousCoordinate _ps = new ContinuousCoordinate();
	
	// random
	private ContinuousCoordinate _pr = new ContinuousCoordinate();
	
	// the set of biomass particle neighbors
	private BiomassParticleSet _neighbors = new BiomassParticleSet();
	
	// ContinuousCoordinates of neighboring voxels
	private ContinuousCoordinate _posX = new ContinuousCoordinate();
	
	private ContinuousCoordinate _negX = new ContinuousCoordinate();
	
	private ContinuousCoordinate _posY = new ContinuousCoordinate();
	
	private ContinuousCoordinate _negY = new ContinuousCoordinate();
	
	private ContinuousCoordinate _posZ = new ContinuousCoordinate();
	
	private ContinuousCoordinate _negZ = new ContinuousCoordinate();
	
	// chemoattractant concentration at cell center
	private float _cc;
	
	// fractional chemoattractant gradient across cell body
	private float _fg;
	
	// chi factors
	private float _chiG;
	
	private float _chiR;
	
	// ContinuousCoordinate of gradient
	private ContinuousCoordinate _gradient = new ContinuousCoordinate();
	
	/**
	 * @param s	the biomass species for a zonal Vicsek SPP
	 */
	public ZonalVicsekBiomassParticle(ZonalVicsekBiomassSpecies s) {
		super(s);
		_species = s;
		_polarization = ContinuousCoordinateOps.getRandomUnitVector();
		_stopOnPeriodicBoundaryCrossing = 
			s._stopOnPeriodicBoundaryCrossing;
	}
	
	/**
	 * sum the gradient contributions to get the polarization this turn
	 */
	private void updatePolarization() {
		// update the polarization terms
		updatePg();
		updatePs();
		updatePr();
		// update chi factors
		updateChis();
		_polarization.reset();
		_polarization.add(ContinuousCoordinateOps.newScale(_pg, _species._wg * _chiG));
		_polarization.add(ContinuousCoordinateOps.newScale(_ps, _species._ws));
		_polarization.add(ContinuousCoordinateOps.newScale(_pr, _species._wr * _chiR));
		// normalize the polarization vector
		if (_polarization.x * _polarization.x 
				+ _polarization.y * _polarization.y 
				+ _polarization.z * _polarization.z > 0) {
			_polarization.normalize();
		}
	}
	
	private void updatePg() {
		_pg.reset();
		_pg.add(_gradient);
		if (_pg.x * _pg.x + _pg.y * _pg.y + _pg.z * _pg.z > 0) {
			_pg.normalize();
		}
		_pg.scale(-1f);
	}
	
	/**
	 * sum polarization contribution from all neighbors and update
	 */
	private void updatePs() {
		_ps.reset();
		while (_neighbors.hasNext()) {
			ZonalVicsekBiomassParticle neib = 
				(ZonalVicsekBiomassParticle) _neighbors.next();
			ContinuousCoordinate dpji = computeNeighborPolarizationContribution(neib);
			_ps.add(dpji);
			if (_ps.x * _ps.x + _ps.y * _ps.y + _ps.z * _ps.z > 0) {
				_ps.normalize();
			}
		}
	}
	
	/**
	 * generate a random polarization vector
	 */
	private void updatePr() {
		_pr = ContinuousCoordinateOps.getRandomUnitVector();
	}
	
	/**
	 * get the polarization contribution from a single neighbor
	 * 
	 * @param neib	a neighboring PolarizableBiomassParticle
	 * @return	the contribution of this neighboring cell
	 */
	private ContinuousCoordinate 
		computeNeighborPolarizationContribution(ZonalVicsekBiomassParticle neib) {
		// calculate the unit direction vector from cell i to cell j
//		ContinuousCoordinate rj = neib.getCenter();
//		ContinuousCoordinate rji = new ContinuousCoordinate(rj.x - _centerX,
//				rj.y - _centerY, rj.z - _centerZ);
//		rji.normalize();
		ContinuousCoordinate pj = neib.getPolarization();
//		float lp = ((pj.x - _polarization.x) * rji.x) + 
//			((pj.y - _polarization.y) * rji.y) +
//			((pj.z - _polarization.z) * rji.z);
//		rji.scale(lp);
//		return rji;
		return pj;
	}
	
	/**
	 * update the chi factors
	 */
	private void updateChis() {
		float hc = _species._co;
		float func = (float) (Math.pow(_cc, hc) / 
				(Math.pow(_species._kc, hc) + Math.pow(_cc, hc)));
		float fung = (float) (Math.pow(_fg, hc) / 
				(Math.pow(_species._kg, hc) + Math.pow(_fg, hc)));
		_chiG = func * fung;
		_chiR = func * (1 - fung);
	}
	
	/**
	 * get the chemoattractant concentration at cell center
	 * TODO This is the production run method.
	 */
	private void updateConcentrationAtCenter() {
		_cc = _species._chemoattractant.getValueAt(this.getCenter());
	}
	
	/**
	 * Update the set of neighboring cells.
	 * (Only homophilic interactions for now.)
	 */
	private void updateNeighbors() {
		_neighbors.reset();
		BiomassParticleSet potentialNeighbors = 
			_particleContainer.getPotentialShovers(this);
		while (potentialNeighbors.hasNext()) {
			BiomassParticle neib = potentialNeighbors.next();
			if (neib.isOfSpecies(_species)) {
				float range = 
					((ZonalVicsekBiomassParticle) neib).getRange() 
					+ this.getRange();
				float distance = ExtraMath.pointDistance(neib.getCenter(), 
						this.getCenter());
				if (range > distance) {
					_neighbors.add(neib);
				}
			}
		}
	}
	
	/**
	 * Update the ContinuousCoordinates of neighboring voxels.
	 * Use periodic boundary conditions in all dimensions.
	 */
	private void updateNeighboringVoxels() {
		float voxelSide = Model.model().getGridElementSide();
		float posX = (_centerX + voxelSide < Model.model().systemSize.x) ?
				_centerX + voxelSide : voxelSide;
		float negX = (_centerX - voxelSide > 0) ? 
				_centerX - voxelSide : Model.model().systemSize.x - voxelSide;
		float posY = (_centerY + voxelSide < Model.model().systemSize.y) ?
				_centerY + voxelSide : voxelSide;
		float negY = (_centerY - voxelSide > 0) ?
				_centerY - voxelSide : Model.model().systemSize.y - voxelSide;
		float posZ = (_centerZ + voxelSide < Model.model().systemSize.z) ?
				_centerZ + voxelSide : voxelSide;
		float negZ = (_centerZ - voxelSide > 0) ?
				_centerZ - voxelSide : Model.model().systemSize.z - voxelSide;
		_posX.x = posX;	_posX.y = _centerY;	_posX.z = _centerZ;
		_negX.x = negX;	_negX.y = _centerY;	_negX.z = _centerZ;
		_posY.x = _centerX;	_posY.y = posY;	_posY.z = _centerZ;
		_negY.x = _centerX;	_negY.y = negY;	_negY.z = _centerZ;
		_posZ.x = _centerX; _posZ.y = _centerY; _posZ.z = posZ;
		_negZ.x = _centerX; _negZ.y = _centerY; _negZ.z = negZ;
	}
	
	/**
	 * Compute the chemoattractant gradient vector across the cell diameter.
	 * Use periodic boundary conditions in all dimensions.
	 * 
	 * TODO This is the production run method.
	 */
	private void updateGradient() {
		float voxelSide = Model.model().getGridElementSide();
		float gradX = (_species._chemoattractant.getValueAt(_posX) - 
				_species._chemoattractant.getValueAt(_negX)) / (2 * voxelSide);
		float gradY = (_species._chemoattractant.getValueAt(_posY) - 
				_species._chemoattractant.getValueAt(_negY)) / (2 * voxelSide);
		float gradZ = (_species._chemoattractant.getValueAt(_posZ) -
				_species._chemoattractant.getValueAt(_negZ)) / (2 * voxelSide);
		_gradient.x = gradX;
		_gradient.y = gradY;
		_gradient.z = gradZ;
		// finally update the fractional gradient
		_fg = (float) Math.sqrt(Math.pow(_gradient.x, 2) 
				+ Math.pow(_gradient.y, 2) 
				+ Math.pow(_gradient.z, 2))
				* 2 * _r / _cc;
	}
	
	/**
	 * Set position of center. 
	 * Overriden for periodic boundary conditions in x dimension.
	 * Overriden to set simulation to finish when a cell crosses periodic
	 * boundary.
	 * 
	 * @param x
	 *            height [micron]
	 * @param y
	 *            y [micron]
	 * @param z
	 *            z [micron]
	 */
	@Override
	public void setCenter(float x, float y, float z) {
		// implement periodic boundaries for x dimension
		float Lxx = _model.systemSize.x;
		if (x >= Lxx) {
			_centerX = x - Lxx;
			if (_stopOnPeriodicBoundaryCrossing)
				Model.model().setFinishSimulation();
		}
		else if(x < 0) {
			_centerX = Lxx + x;
			if (_stopOnPeriodicBoundaryCrossing)
				Model.model().setFinishSimulation();
		}
		else
			_centerX = x;
		// Implement periodic boundaries (side walls)
		float Lyy = _model.systemSize.y;
		if (y >= Lyy) {
			_centerY = y - Lyy;
			if (_stopOnPeriodicBoundaryCrossing)
				Model.model().setFinishSimulation();
		}
		else if (y < 0) {
			_centerY = Lyy + y;
			if (_stopOnPeriodicBoundaryCrossing)
				Model.model().setFinishSimulation();
		}
		else
			_centerY = y;
		float Lzz = _model.systemSize.z;
		if (Lzz == 0)
			// 2D case
			_centerZ = 0;
		else if (z >= Lzz) {
			_centerZ = z - Lzz;
			if (_stopOnPeriodicBoundaryCrossing)
				Model.model().setFinishSimulation();
		}
		else if (z < 0) {
			_centerZ = Lzz + z;
			if (_stopOnPeriodicBoundaryCrossing)
				Model.model().setFinishSimulation();
		}
		else
			_centerZ = z;
		setShovingMapPositionAndAddToShovingMap();
	}
	
	/**
	 * Move the cell by actin remodeling if there are neighbor cells or
	 * chemoattractant present
	 * 
	 * @param t	time step length
	 */
	private void migrate(float t) {
//		float speed = (float) Math.sqrt(_species._motility * t);
		float speed = _species._motility;
		ContinuousCoordinate displacement = 
			ContinuousCoordinateOps.newScale(_polarization, speed * t);
		move(displacement.x, displacement.y, displacement.z);
	}
	
	/**
	 * Decide whether to migrate instead of grow this turn
	 */
	private void makeMigrationDecision() {
		updateConcentrationAtCenter();
		// migration OR growth in a single turn! 
		float func = (float) (Math.pow(_cc, _species._co) / 
				(Math.pow(_species._kc, _species._co) + Math.pow(_cc, _species._co)));
		// throw random number to choose whether to migrate instead of grow
		_migratingThisIteration = (func >= Model.model().getRandom());
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
		float timeStep;
		makeMigrationDecision();
		// if migrating calculate time step based on motility
		if (_migratingThisIteration) {
			timeStep = _r / _species._motility; 
		}
		// otherwise calculate time step based o growth
		else {
			timeStep = super.getMaximumTimeStep();
		}
		// force the model max time step if less than cell-specific time step
//		timeStep = Math.min(Model.model().getMaximumTimeStep(), timeStep);
		return timeStep;
	}
	
	@Override
	public float grow(float t) {
		float m;		
		if (_migratingThisIteration) {
			// update environmental cues
			updateNeighboringVoxels();
			updateNeighbors();
			updateGradient();
			updateChis();
			// update polarization vector
			updatePolarization();
			// TODO code to account for shoving in polarization
			migrate(t);
			/*
			 * Here we assume that in a migrating cell all biomass reactions
			 * are stopped (growth, decay etc.). This probably will need to be
			 * changed.
			 */
		}
		else {
			_polarization.reset();
		}
		m = super.grow(t);
		return m;
	}
	
	@Override
	/**
	 * Override the superclass method to draw an orange capsule if cell is
	 * moving
	 */
	public Color getColorCapsule() {
		if (_migratingThisIteration) {
			return Color.ORANGE;
		}
		else {
			return super.getColorCapsule();
		}
	}
	
	/**
	 * Override toString() to additionally write the polarization vector
	 * components to the output file.
	 * 
	 * @see java.lang.Object#toString()
	 */
	@Override
	public String toString() {
		Color colorCore = getColorCore();
		Color colorCapsule = getColorCapsule();
		return _centerX + "\t" + _centerY + "\t" + _centerZ + "\t"
				+ getCoreRadius() + "\t" + colorCore.getRed() + "\t"
				+ colorCore.getGreen() + "\t" + colorCore.getBlue() + "\t" + _r
				+ "\t" + colorCapsule.getRed() + "\t" + colorCapsule.getGreen()
				+ "\t" + colorCapsule.getBlue() 
				+ "\t" + _polarization.x + "\t" + _polarization.y + "\t"
				+ _polarization.z + "\t" + _species.getName() + "\n";
	}
	
	/**
	 * public method for getting the adhesive radius
	 * 
	 * @return	the adhesive radius
	 */
	public float getRange() {
		return _r * _species._ra;
	}
	
	/**
	 * public method for getting the polarization
	 * 
	 * @return	the polarization
	 */
	public ContinuousCoordinate getPolarization() {
		return _polarization;
	}

}

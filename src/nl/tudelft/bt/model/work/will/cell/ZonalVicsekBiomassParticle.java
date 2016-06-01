/**
 * 
 */
package nl.tudelft.bt.model.work.will.cell;

import nl.tudelft.bt.model.ContinuousCoordinate;
import nl.tudelft.bt.model.particlebased.BiomassParticleSet;
import nl.tudelft.bt.model.util.ContinuousCoordinateOps;
import nl.tudelft.bt.model.work.will.MutatorBiomassParticle;
import nl.tudelft.bt.model.work.will.MutatorBiomassSpecies;
import nl.tudelft.bt.model.work.will.PolarizableBiomassParticle;

/**
 * @author wkc
 *
 */
public class ZonalVicsekBiomassParticle extends MutatorBiomassParticle {
	
	private ZonalVicsekBiomassSpecies _species;

	// stop simulation on periodic boundary crossing?
	private boolean _stopOnPeriodicBoundaryCrossing;
	
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
	private void updateTotalPolarization() {
		// update the polarization terms
		updateGradientSensingPolarization();
		updateSocialPolarization();
		updateRandomWalkPolarization();
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
	
	private void updateGradientSensingPolarization() {
		_pg.reset();
		_pg.add(_gradient);
		if (_pg.x * _pg.x + _pg.y * _pg.y + _pg.z * _pg.z > 0) {
			_pg.normalize();
		}
	}
	
	/**
	 * sum polarization contribution from all neighbors and update
	 */
	private void updateSocialPolarization() {
		_ps.reset();
		while (_neighbors.hasNext()) {
			PolarizableBiomassParticle neib = 
				(PolarizableBiomassParticle) _neighbors.next();
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
	private void updateRandomWalkPolarization() {
		_pr = ContinuousCoordinateOps.getRandomUnitVector();
	}
	
	/**
	 * get the polarization contribution from a single neighbor
	 * 
	 * @param neib	a neighboring PolarizableBiomassParticle
	 * @return	the contribution of this neighboring cell
	 */
	private ContinuousCoordinate 
		computeNeighborPolarizationContribution(PolarizableBiomassParticle neib) {
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

}

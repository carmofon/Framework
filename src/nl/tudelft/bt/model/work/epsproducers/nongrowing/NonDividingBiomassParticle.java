package nl.tudelft.bt.model.work.epsproducers.nongrowing;

import nl.tudelft.bt.model.particlebased.BiomassParticle;

/**
 * Implements surface motility into a biomass particle
 * 
 * @author Joao Xavier (jxavier@cgr.harvard.edu) - Oct 12, 2006
 */
public class NonDividingBiomassParticle extends BiomassParticle {
	private boolean _isEpsOnly = false;

	public NonDividingBiomassParticle(NonDividingBiomassSpecies s) {
		super(s);
	}

	@Override
	public boolean willDivide() {
		// divide only if it is EPS particle
		if (_isEpsOnly)
			return super.willDivide();
		// otherwise returs false
		return false;
	}

	
	
	@Override
	public BiomassParticle excreteEps() {
		NonDividingBiomassParticle p = (NonDividingBiomassParticle)super.excreteEps();
		p._isEpsOnly = true;
		return p;
	}

	@Override
	protected Object clone() throws CloneNotSupportedException {
		// TODO Auto-generated method stub
		return super.clone();
	}
}
